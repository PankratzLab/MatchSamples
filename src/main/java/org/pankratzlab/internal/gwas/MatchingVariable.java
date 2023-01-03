package org.pankratzlab.internal.gwas;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.pankratzlab.common.ext;

public class MatchingVariable {
  final String prettyName;
  final String headerName;
  final boolean isBinary;

  private int headerIndex = -1;

  // for continuous variables
  private int casesSeen = 0;
  private int controlsSeen = 0;

  private double caseSum = 0;
  private double controlSum = 0;

  private double caseAvg;
  private double controlAvg;

  // for binary variables
  Map<String, Integer> caseValues = new HashMap<>();
  private Map<String, List<Integer>> previousControlValues = new HashMap<>();

  private int totalSeen = 0;
  private int matchesSeen = 0;

  private double concordance;

  public MatchingVariable(String headerName, boolean isBinary) {
    this.headerName = headerName;
    this.prettyName = headerName;
    this.isBinary = isBinary;
  }

  public int findIndexInHeader(String[] header) {
    headerIndex = ext.indexOfStr(headerName, header);
    return headerIndex;
  }

  public int getHeaderIndex() {
    return headerIndex;
  }

  public void addToCaseSum(double value) {
    caseSum += value;
    casesSeen += 1;
  }

  public void addToControlSum(double value) {
    controlSum += value;
    controlsSeen += 1;
  }

  public void setCaseValue(String caseId, int value) {
    caseValues.put(caseId, value);
    if (previousControlValues.containsKey(caseId)) {
      for (int v : previousControlValues.get(caseId)) {
        if (v == value) {
          matchesSeen += 1;
        }
        totalSeen += 1;
      }
    }
    previousControlValues.remove(caseId);
  }

  public void recordControlValue(String caseId, int value) {
    if (caseValues.containsKey(caseId)) {
      if (value == caseValues.get(caseId)) {
        matchesSeen += 1;
      }
      totalSeen += 1;
    } else {
      if (previousControlValues.containsKey(caseId)) {
        previousControlValues.get(caseId).add(value);
      } else {
        previousControlValues.put(caseId, new ArrayList<>(List.of(value)));
      }
    }
  }

  public double getCaseAvg() {
    this.caseAvg = caseSum / casesSeen;
    return caseAvg;
  }

  public double getControlAvg() {
    controlAvg = controlSum / controlsSeen;
    return controlAvg;
  }

  public double getConcordance() {
    concordance = (double) matchesSeen / (double) totalSeen;
    return concordance;
  }
}
