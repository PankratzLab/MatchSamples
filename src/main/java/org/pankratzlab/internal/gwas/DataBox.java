package org.pankratzlab.internal.gwas;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import org.pankratzlab.common.ext;
import org.pankratzlab.common.stats.LogisticRegression;
import org.pankratzlab.common.stats.RegressionModel;

public class DataBox {

  private final MatchingVariable[] matchingVariables;

  // keep a map of MV -> index in array for O(1) lookup later
  private final Map<MatchingVariable, Integer> matchingVariableIndexMap;

  private final int totalSampleCount;
  private final double[][] data;

  private final Map<String, String> controlCasePairings;
  private final Map<String, Integer> sampleIndexById;
  private final String[] sampleIdByIndex;
  private final Set<String> caseIds;
  private final Set<String> controlIds;

  private int nextSampleIndex = 0;

  private double[] concordances;
  private boolean concordancesComputed = false;

  private double[] caseAverages;
  private double[] controlAverages;
  private boolean averagesComputed = false;

  private double[] univariatePValues;
  private boolean univariatePValuesComputed = false;

  private double[] multivariatePValues;
  private boolean multivariatePValuesComputed = false;

  private final Logger log = Logger.getAnonymousLogger();

  public DataBox(MatchingVariable[] matchingVariables, Map<String, String> controlCasePairings) {
    this.matchingVariables = matchingVariables;

    this.matchingVariableIndexMap = new HashMap<>();
    for (int i = 0; i < this.matchingVariables.length; i++) {
      MatchingVariable mv = this.matchingVariables[i];
      matchingVariableIndexMap.put(mv, i);
      mv.setDataBox(this);
    }

    this.controlCasePairings = controlCasePairings;
    this.caseIds = new HashSet<>(controlCasePairings.values());
    this.controlIds = controlCasePairings.keySet();
    this.totalSampleCount = caseIds.size() + controlIds.size();

    this.data = new double[totalSampleCount][matchingVariables.length];
    this.sampleIndexById = new HashMap<>(totalSampleCount);
    this.sampleIdByIndex = new String[totalSampleCount];
  }

  public void recordData(String[] line) {
    String sampleId = line[0];
    if (caseIds.contains(sampleId) || controlIds.contains(sampleId)) {
      if (sampleIndexById.containsKey(sampleId)) {
        throw new IllegalStateException("Found the same sample ID twice: " + sampleId);
      }
      sampleIndexById.put(sampleId, nextSampleIndex);
      sampleIdByIndex[nextSampleIndex] = sampleId;

      for (int i = 0; i < matchingVariables.length; i++) {
        MatchingVariable mv = matchingVariables[i];
        double value = Double.parseDouble(line[mv.getHeaderIndex()]);
        mv.checkBinary(value);
        data[nextSampleIndex][i] = value;
      }
      nextSampleIndex++;
    }
  }

  public void computeConcordances() {
    int[] matchCounts = new int[matchingVariables.length];
    this.concordances = new double[matchingVariables.length];

    // iterate over samples
    for (int si = 0; si < totalSampleCount; si++) {
      if (caseIds.contains(sampleIdByIndex[si])) {
        // we want to look at controls only
        continue;
      }
      // From here onwards we know the sample in question is a control
      int caseIndex = getCaseIndexForControlIndex(si);
      // iterate over variables
      for (int vi = 0; vi < matchCounts.length; vi++) {
        if (matchingVariables[vi].isBinary()) {
          double controlValue = data[si][vi];
          double caseValue = data[caseIndex][vi];
          if (controlValue == caseValue) {
            matchCounts[vi]++;
          }
        }
      }
    }

    for (int i = 0; i < concordances.length; i++) {
      if (matchingVariables[i].isBinary()) {
        // divide by the number of controls
        concordances[i] = (double) matchCounts[i] / controlCasePairings.size();
      } else {
        concordances[i] = Double.NaN;
      }
    }
    concordancesComputed = true;
  }

  public void computeAverages() {
    this.caseAverages = new double[matchingVariables.length];
    this.controlAverages = new double[matchingVariables.length];

    // iterate over samples
    for (int si = 0; si < totalSampleCount; si++) {
      // iterate over variables
      for (int vi = 0; vi < matchingVariables.length; vi++) {
        if (matchingVariables[vi].isContinuous()) {
          if (caseIds.contains(sampleIdByIndex[si])) {
            caseAverages[vi] += data[si][vi];
          } else {
            controlAverages[vi] += data[si][vi];
          }
        }
      }
    }
    for (int i = 0; i < matchingVariables.length; i++) {
      if (matchingVariables[i].isContinuous()) {
        // divide by number of cases and controls, respectively
        caseAverages[i] /= caseIds.size();
        controlAverages[i] /= controlCasePairings.size();
      } else {
        caseAverages[i] = controlAverages[i] = Double.NaN;
      }
    }
    averagesComputed = true;
  }

  public void computeUnivariateP() {
    double[] deps = Arrays.stream(sampleIdByIndex)
                          .mapToDouble(sampleId -> caseIds.contains(sampleId) ? 1 : 0).toArray();

    this.univariatePValues = new double[matchingVariables.length];

    for (int mvIndex = 0; mvIndex < matchingVariables.length; mvIndex++) {
      MatchingVariable mv = matchingVariables[mvIndex];

      double[] indeps = new double[totalSampleCount];

      for (int si = 0; si < totalSampleCount; si++) {
        indeps[si] = data[si][mvIndex];
      }

      RegressionModel model = new LogisticRegression(deps, indeps);
      univariatePValues[mvIndex] = model.getOverallSig();
    }
    univariatePValuesComputed = true;
  }

  public void computeMultivariateP() {
    this.multivariatePValues = new double[matchingVariables.length];
    // dependent variables: we just have one, case/control status
    // this is represented as an array of 1s for cases and 0s for controls
    double[] deps = Arrays.stream(sampleIdByIndex)
                          .mapToDouble(sampleId -> caseIds.contains(sampleId) ? 1 : 0).toArray();

    String[] indepVariableNames = Arrays.stream(matchingVariables)
                                        .map(matchingVariable -> matchingVariable.headerName)
                                        .toArray(String[]::new);

    RegressionModel model = new LogisticRegression(deps, data, indepVariableNames, false, true);

    for (MatchingVariable mv : matchingVariables) {
      int indexInModelSigs = ext.indexOfStr(mv.headerName, model.getVarNames());
      double p = Double.NaN;
      if (indexInModelSigs == -1) {
        log.warning("No multivariate p calculated for " + mv.headerName
                    + ". Value will be set to NaN.");
      } else {
        p = model.getSigs()[indexInModelSigs];
      }
      multivariatePValues[matchingVariableIndexMap.get(mv)] = p;

      if (Double.isNaN(p)) {
        log.warning("Found multivariate p value for " + mv.headerName
                    + " is NaN. This may be caused by colinearity or other malformations in the phenotype data.");
      }
    }
    multivariatePValuesComputed = true;
  }

  public void computeConcordancesIfNecessary() {
    if (!concordancesComputed) {
      computeConcordances();
    }
  }

  public void computeAveragesIfNecessary() {
    if (!averagesComputed) {
      computeAverages();
    }
  }

  public void computeUnivariatePIfNecessary() {
    if (!univariatePValuesComputed) {
      computeUnivariateP();
    }
  }

  public void computeMultivariatePIfNecessary() {
    if (!multivariatePValuesComputed) {
      computeMultivariateP();
    }
  }

  public double getConcordance(MatchingVariable matchingVariable) {
    matchingVariable.mustBeBinary();
    computeConcordancesIfNecessary();
    int mvIndex = matchingVariableIndexMap.get(matchingVariable);
    return concordances[mvIndex];
  }

  public double getCaseAvg(MatchingVariable matchingVariable) {
    matchingVariable.cantBeBinary();
    computeAveragesIfNecessary();
    int mvIndex = matchingVariableIndexMap.get(matchingVariable);
    return caseAverages[mvIndex];
  }

  public double getControlAvg(MatchingVariable matchingVariable) {
    matchingVariable.cantBeBinary();
    computeAveragesIfNecessary();
    int mvIndex = matchingVariableIndexMap.get(matchingVariable);
    return controlAverages[mvIndex];
  }

  public double getUnivariateP(MatchingVariable matchingVariable) {
    computeUnivariatePIfNecessary();
    int combinedMvIndex = matchingVariableIndexMap.get(matchingVariable);
    return univariatePValues[combinedMvIndex];
  }

  public double getMultivariateP(MatchingVariable matchingVariable) {
    computeMultivariatePIfNecessary();
    int combinedMvIndex = matchingVariableIndexMap.get(matchingVariable);
    return multivariatePValues[combinedMvIndex];
  }

  private int getCaseIndexForControlIndex(int controlIndex) {
    String controlId = sampleIdByIndex[controlIndex];
    String caseId = controlCasePairings.get(controlId);
    if (sampleIndexById.containsKey(caseId)) {
      return sampleIndexById.get(caseId);
    } else {
      throw new IllegalStateException("No index exists for this case: " + caseId
                                      + ", which is the control for " + controlId);
    }
  }
}
