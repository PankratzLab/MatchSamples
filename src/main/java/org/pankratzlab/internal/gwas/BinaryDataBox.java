package org.pankratzlab.internal.gwas;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class BinaryDataBox {
  private final Map<String, int[]> valuesBySampleId;
  private final Map<MatchingVariable, Integer> indexByMatchingVariable;
  private final Map<String, String> controlCasePairings;
  private final Set<String> caseIds;

  private double[] concordances;
  private boolean computed = false;
  private final int matchingVariableCount;
  private final int expectedCaseValuesToSet;
  private int actualCaseValuesSet = 0;

  public BinaryDataBox(MatchingVariable[] matchingVariables, int totalSamples,
                       Map<String, String> controlCasePairings, Set<String> caseIds) {
    this.matchingVariableCount = matchingVariables.length;
    this.indexByMatchingVariable = new HashMap<>(matchingVariables.length);
    for (int i = 0; i < matchingVariables.length; i++) {
      MatchingVariable mv = matchingVariables[i];
      if (!mv.isBinary) {
        throw new IllegalArgumentException("Cannot create a binary data box with a non-binary matching variable");
      }
      indexByMatchingVariable.put(mv, i);
      mv.setBinaryDataBox(this);
    }

    this.valuesBySampleId = new HashMap<>(totalSamples);
    this.controlCasePairings = controlCasePairings;
    this.caseIds = caseIds;

    expectedCaseValuesToSet = caseIds.size() * matchingVariableCount;
  }

  public void recordValue(String sampleId, MatchingVariable matchingVariable, int value) {
    int mvIndex = indexByMatchingVariable.get(matchingVariable);
    if (valuesBySampleId.containsKey(sampleId)) {
      valuesBySampleId.get(sampleId)[mvIndex] = value;
    } else {
      int[] values = new int[matchingVariableCount];
      values[mvIndex] = value;
      valuesBySampleId.put(sampleId, values);
    }
    if (caseIds.contains(sampleId)) {
      actualCaseValuesSet++;
    }
  }

  public void computeConcordances() {
    if (actualCaseValuesSet < expectedCaseValuesToSet) {
      throw new IllegalStateException("Expected to set " + expectedCaseValuesToSet
                                      + " case values, but only " + actualCaseValuesSet
                                      + " have been set.");
    }
    int[] matchCounts = new int[matchingVariableCount];
    for (String controlId : controlCasePairings.keySet()) {
      String caseId = controlCasePairings.get(controlId);
      for (int i = 0; i < matchingVariableCount; i++) {
        if (valuesBySampleId.get(controlId)[i] == valuesBySampleId.get(caseId)[i]) {
          matchCounts[i]++;
        }
      }
    }
    concordances = new double[matchingVariableCount];
    for (int i = 0; i < concordances.length; i++) {
      // we want the percentage of controls that match their cases, so only divide
      // by the number of controls
      concordances[i] = (double) matchCounts[i] / controlCasePairings.size();
    }
    computed = true;
  }

  public double getConcordance(MatchingVariable matchingVariable) {
    if (!computed) {
      computeConcordances();
    }
    int mvIndex = indexByMatchingVariable.get(matchingVariable);
    return this.concordances[mvIndex];
  }
}
