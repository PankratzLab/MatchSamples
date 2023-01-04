package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.pankratzlab.common.Files;
import org.pankratzlab.common.PSF;
import org.pankratzlab.kdmatch.KDMatch;

public class REval {
  final MatchingVariable[] matchingVariables;
  final File statusFile;
  final File phenotypeFile;
  final Map<String, String> controlCasePairings;
  final Set<String> cases;
  final Set<String> controls;

  public REval(MatchingVariable[] matchingVariables, File statusFile,
               File phenotypeFile) throws IOException {
    this.matchingVariables = matchingVariables;

    if (!statusFile.exists()) {
      throw new IllegalArgumentException("Provided status file does not exist.");
    }
    if (!phenotypeFile.exists()) {
      throw new IllegalArgumentException("Provided phenotype file does not exist.");
    }

    this.statusFile = statusFile;
    this.phenotypeFile = phenotypeFile;

    this.controlCasePairings = readPairings();
    this.controls = controlCasePairings.keySet();
    this.cases = new HashSet<>(controlCasePairings.values());
  }

  public void evaluate() throws IOException {
    BufferedReader phenoReader = Files.getAppropriateReader(phenotypeFile.toString());
    String[] phenoHeader = splitTsvLine(phenoReader.readLine());

    // figure out which column corresponds to each matching variable
    Map<Integer, MatchingVariable> matchingVarIndexMap = new HashMap<>();
    for (MatchingVariable mv : matchingVariables) {
      int i = mv.findIndexInHeader(phenoHeader);
      if (i == -1) {
        throw new IllegalStateException("Unable to find column for matching variable; expected header name: "
                                        + mv.headerName);
      }
      MatchingVariable previousMV = matchingVarIndexMap.putIfAbsent(i, mv);
      if (previousMV != null) {
        System.out.println("WARNING - Found multiple matching variables for the same header name: "
                           + previousMV.headerName);
      }
    }

    String[] line;
    String inputLine = phenoReader.readLine();
    while (inputLine != null) {
      line = splitTsvLine(inputLine);
      String sampleId = line[0];
      if (cases.contains(sampleId) || controls.contains(sampleId)) {
        for (int matchVarIndex : matchingVarIndexMap.keySet()) {
          MatchingVariable mv = matchingVarIndexMap.get(matchVarIndex);
          if (mv.isBinary) {
            if (cases.contains(sampleId)) {
              mv.setCaseValue(sampleId, Integer.parseInt(line[matchVarIndex]));
            } else {
              mv.recordControlValue(controlCasePairings.get(sampleId),
                                    Integer.parseInt(line[matchVarIndex]));
            }
          } else {
            double value = Double.parseDouble(line[matchVarIndex]);
            if (cases.contains(sampleId)) {
              mv.addToCaseSum(value);
            } else {
              mv.addToControlSum(value);
            }
          }
        }
      }
      inputLine = phenoReader.readLine();
    }
  }

  private Map<String, String> readPairings() throws IOException {
    // id -> matched case
    Map<String, String> pairings = new HashMap<>();
    String expectedHeader = KDMatch.STATUS_FILE_HEADER;
    BufferedReader reader = Files.getAppropriateReader(statusFile.toString());
    String actualHeader = reader.readLine().strip();
    if (!actualHeader.equals(expectedHeader)) {
      throw new IllegalStateException("Status file header does not match expected header.");
    }

    String[] line;
    String inputLine = reader.readLine();
    while (inputLine != null) {
      line = splitTsvLine(inputLine);
      // We don't want any (case -> itself) pairs
      // we don't care about a case's relationship to itself
      boolean isNotCase = Integer.parseInt(line[1]) != 1;
      if (isNotCase) {
        pairings.put(line[0], line[2]);
      }
      inputLine = reader.readLine();
    }
    return pairings;
  }

  private static String[] splitTsvLine(String inputLine) {
    return inputLine.trim().split(PSF.Regex.GREEDY_WHITESPACE);
  }
}
