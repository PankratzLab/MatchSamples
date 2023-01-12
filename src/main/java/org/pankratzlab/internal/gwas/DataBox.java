package org.pankratzlab.internal.gwas;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.pankratzlab.common.ext;
import org.pankratzlab.common.stats.LogisticRegression;
import org.pankratzlab.common.stats.RegressionModel;

public class DataBox {

  private final MatchingVariable[] continuousMatchingVariables;
  private final MatchingVariable[] binaryMatchingVariables;

  // keep a map of MV -> index in array for O(1) lookup later
  // This stores each MV's index AMONG the other MVs OF ITS TYPE.
  private final Map<MatchingVariable, Integer> matchingVariableIndexMap;
  // This is an index to be used when we need a unique index of each MV among ALL the MVs.
  // The structure is simple: continuous MVs in their original order, then binary MVs in their
  // original order.
  private final Map<MatchingVariable, Integer> combinedVariableIndexMap;

  private final int totalSampleCount;
  private final int[][] binaryData;
  private final double[][] continuousData;

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

  public DataBox(MatchingVariable[] allMatchingVariables, Map<String, String> controlCasePairings) {
    this(Arrays.stream(allMatchingVariables).filter(mv -> !mv.isBinary)
               .toArray(MatchingVariable[]::new),
         Arrays.stream(allMatchingVariables).filter(mv -> mv.isBinary)
               .toArray(MatchingVariable[]::new),
         controlCasePairings);
  }

  public DataBox(MatchingVariable[] continuousMatchingVariables,
                 MatchingVariable[] binaryMatchingVariables,
                 Map<String, String> controlCasePairings) {

    this.continuousMatchingVariables = continuousMatchingVariables;
    this.binaryMatchingVariables = binaryMatchingVariables;

    // Initialize MV index maps and set data box fields in each MV.
    this.matchingVariableIndexMap = new HashMap<>();
    this.combinedVariableIndexMap = new HashMap<>();
    for (int i = 0; i < continuousMatchingVariables.length; i++) {
      MatchingVariable mv = continuousMatchingVariables[i];
      matchingVariableIndexMap.put(mv, i);
      combinedVariableIndexMap.put(mv, i);
      mv.setDataBox(this);
    }
    for (int i = 0; i < binaryMatchingVariables.length; i++) {
      MatchingVariable mv = binaryMatchingVariables[i];
      matchingVariableIndexMap.put(mv, i);
      combinedVariableIndexMap.put(mv, continuousMatchingVariables.length + i);
      mv.setDataBox(this);
    }

    this.controlCasePairings = controlCasePairings;
    this.caseIds = new HashSet<>(controlCasePairings.values());
    this.controlIds = controlCasePairings.keySet();
    this.totalSampleCount = caseIds.size() + controlIds.size();

    this.binaryData = new int[totalSampleCount][binaryMatchingVariables.length];
    this.continuousData = new double[totalSampleCount][continuousMatchingVariables.length];

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

      for (int i = 0; i < binaryMatchingVariables.length; i++) {
        MatchingVariable mv = binaryMatchingVariables[i];
        binaryData[nextSampleIndex][i] = Integer.parseInt(line[mv.getHeaderIndex()]);
      }
      for (int i = 0; i < continuousMatchingVariables.length; i++) {
        MatchingVariable mv = continuousMatchingVariables[i];
        continuousData[nextSampleIndex][i] = Double.parseDouble(line[mv.getHeaderIndex()]);
      }
      nextSampleIndex++;
    }
  }

  public void computeConcordances() {
    int[] matchCounts = new int[binaryMatchingVariables.length];
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
        int controlValue = binaryData[si][vi];
        int caseValue = binaryData[caseIndex][vi];
        if (controlValue == caseValue) {
          matchCounts[vi]++;
        }
      }
    }
    this.concordances = new double[binaryMatchingVariables.length];
    for (int i = 0; i < concordances.length; i++) {
      // divide by the number of controls
      concordances[i] = (double) matchCounts[i] / controlCasePairings.size();
    }
    concordancesComputed = true;
  }

  public void computeAverages() {
    this.caseAverages = new double[continuousMatchingVariables.length];
    this.controlAverages = new double[continuousMatchingVariables.length];
    // iterate over samples
    for (int si = 0; si < totalSampleCount; si++) {
      // iterate over variables
      for (int vi = 0; vi < continuousMatchingVariables.length; vi++) {
        if (caseIds.contains(sampleIdByIndex[si])) {
          caseAverages[vi] += continuousData[si][vi];
        } else {
          controlAverages[vi] += continuousData[si][vi];
        }
      }
    }
    for (int i = 0; i < continuousMatchingVariables.length; i++) {
      // divide by number of cases and controls, respectively
      caseAverages[i] /= caseIds.size();
      controlAverages[i] /= controlCasePairings.size();
    }
    averagesComputed = true;
  }

  public void computeUnivariateP() {
    double[] deps = Arrays.stream(sampleIdByIndex)
                          .mapToDouble(sampleId -> caseIds.contains(sampleId) ? 1 : 0).toArray();
    this.univariatePValues = new double[matchingVariableIndexMap.size()];

    for (MatchingVariable mv : matchingVariableIndexMap.keySet()) {
      double[] indeps = new double[totalSampleCount];
      int mvIndex = matchingVariableIndexMap.get(mv);
      if (mv.isBinary) {
        for (int si = 0; si < totalSampleCount; si++) {
          indeps[si] = binaryData[si][mvIndex];
        }
      } else {
        for (int si = 0; si < totalSampleCount; si++) {
          indeps[si] = continuousData[si][mvIndex];
        }
      }
      RegressionModel model = new LogisticRegression(deps, indeps);
      univariatePValues[combinedVariableIndexMap.get(mv)] = model.getOverallSig();
    }
    univariatePValuesComputed = true;
  }

  public void computeMultivariateP() {
    this.multivariatePValues = new double[combinedVariableIndexMap.size()];
    // dependant variables: we just have one, case/control status
    // this is represented as an array of 1s for cases and 0s for controls
    double[] deps = Arrays.stream(sampleIdByIndex)
                          .mapToDouble(sampleId -> caseIds.contains(sampleId) ? 1 : 0).toArray();

    double[][] indeps = new double[totalSampleCount][combinedVariableIndexMap.size()];
    for (int si = 0; si < totalSampleCount; si++) {
      for (MatchingVariable mv : combinedVariableIndexMap.keySet()) {
        if (mv.isBinary) {
          indeps[si][combinedVariableIndexMap.get(mv)] = binaryData[si][matchingVariableIndexMap.get(mv)];
        } else {
          indeps[si][combinedVariableIndexMap.get(mv)] = continuousData[si][matchingVariableIndexMap.get(mv)];
        }
      }
    }

    String[] indepVariableNames = new String[combinedVariableIndexMap.size()];
    for (MatchingVariable mv : combinedVariableIndexMap.keySet()) {
      int combinedIndex = combinedVariableIndexMap.get(mv);
      indepVariableNames[combinedIndex] = mv.headerName;
    }

    RegressionModel model = new LogisticRegression(deps, indeps, indepVariableNames, false, true);

    for (MatchingVariable mv : combinedVariableIndexMap.keySet()) {
      int indexInModelSigs = ext.indexOfStr(mv.headerName, model.getVarNames());
      multivariatePValues[combinedVariableIndexMap.get(mv)] = model.getSigs()[indexInModelSigs];
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
    int combinedMvIndex = combinedVariableIndexMap.get(matchingVariable);
    return univariatePValues[combinedMvIndex];
  }

  public double getMultivariateP(MatchingVariable matchingVariable) {
    computeMultivariatePIfNecessary();
    int combinedMvIndex = combinedVariableIndexMap.get(matchingVariable);
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
