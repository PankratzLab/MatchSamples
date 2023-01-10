package org.pankratzlab.internal.gwas;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class ContinuousDataBox {

  public final Map<MatchingVariable, Integer> indexByMatchingVariable;
  private final int totalSamples;
  private final Set<String> caseIds;
  private Map<String, double[]> sampleData;
  private final int casesCount;
  private final int controlsCount;

  private double[] caseAverages;
  private double[] controlAverages;
  private boolean computed = false;

  public ContinuousDataBox(MatchingVariable[] matchingVariables, int totalSamples,
                           Map<String, String> controlCasePairings, Set<String> caseIds) {
    this.indexByMatchingVariable = new HashMap<>(matchingVariables.length);
    for (int i = 0; i < matchingVariables.length; i++) {
      MatchingVariable mv = matchingVariables[i];
      if (mv.isBinary) {
        throw new IllegalArgumentException("Cannot create a continuous data box with a binary matching variable");
      }
      indexByMatchingVariable.put(mv, i);
      mv.setContinuousDataBox(this);
    }

    this.totalSamples = totalSamples;
    this.sampleData = new HashMap<>(totalSamples);
    this.caseIds = caseIds;
    this.casesCount = caseIds.size();
    this.controlsCount = controlCasePairings.size();
  }

  public void recordValue(MatchingVariable matchingVariable, String sampleId, double value) {
    int mvIndex = indexByMatchingVariable.get(matchingVariable);
    if (sampleData.containsKey(sampleId)) {
      sampleData.get(sampleId)[mvIndex] = value;
    } else {
      double[] data = new double[indexByMatchingVariable.size()];
      data[mvIndex] = value;
      sampleData.put(sampleId, data);
    }
  }

  public void computeAverages() {
    caseAverages = new double[indexByMatchingVariable.size()];
    controlAverages = new double[indexByMatchingVariable.size()];

    // iterate through samples
    for (Map.Entry<String, double[]> entry : sampleData.entrySet()) {
      boolean isCase = caseIds.contains(entry.getKey());
      // iterate through matching variables for each sample
      for (int i = 0; i < indexByMatchingVariable.size(); i++) {
        if (isCase) {
          caseAverages[i] += entry.getValue()[i];
        } else {
          controlAverages[i] += entry.getValue()[i];
        }
      }
    }
    for (int mvi = 0; mvi < indexByMatchingVariable.size(); mvi++) {
      caseAverages[mvi] /= casesCount;
      controlAverages[mvi] /= controlsCount;
    }
    computed = true;
  }

  private void computeAveragesIfNecessary() {
    if (!computed) {
      computeAverages();
    }
  }

  public double getCaseAvg(MatchingVariable matchingVariable) {
    computeAveragesIfNecessary();
    int mvIndex = indexByMatchingVariable.get(matchingVariable);
    return caseAverages[mvIndex];
  }

  public double getControlAvg(MatchingVariable matchingVariable) {
    computeAveragesIfNecessary();
    int mvIndex = indexByMatchingVariable.get(matchingVariable);
    return controlAverages[mvIndex];
  }
}
