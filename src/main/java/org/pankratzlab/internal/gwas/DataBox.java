package org.pankratzlab.internal.gwas;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class DataBox {

  private final MatchingVariable[] continuousMatchingVariables;
  private final MatchingVariable[] binaryMatchingVariables;
  // keep a map of MV -> index in array for O(1) lookup later
  private final Map<MatchingVariable, Integer> matchingVariableIndexMap;

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

    // Initialize MV index map and set data box fields in each MV.
    this.matchingVariableIndexMap = new HashMap<>();
    for (int i = 0; i < continuousMatchingVariables.length; i++) {
      MatchingVariable mv = continuousMatchingVariables[i];
      matchingVariableIndexMap.put(mv, i);
      mv.setDataBox(this);
    }
    for (int i = 0; i < binaryMatchingVariables.length; i++) {
      MatchingVariable mv = binaryMatchingVariables[i];
      matchingVariableIndexMap.put(mv, i);
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
