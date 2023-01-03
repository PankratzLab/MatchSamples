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
    // todo: don't use assert
    assert statusFile.exists();
    assert phenotypeFile.exists();

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
      assert i != -1; // todo: checks that all the matching variables end up in the map
      matchingVarIndexMap.put(i, mv);
      // todo: there should be some check that no two matching variables share an index
    }

    String[] line;
    String inputLine = phenoReader.readLine();
    while (inputLine != null) {
      line = splitTsvLine(inputLine);
      String sampleId = line[0];
      if (cases.contains(sampleId) || controls.contains(sampleId)) {
        for (int i : matchingVarIndexMap.keySet()) {
          MatchingVariable mv = matchingVarIndexMap.get(i);
          if (mv.isBinary) {
            if (cases.contains(sampleId)) {
              mv.setCaseValue(sampleId, Integer.parseInt(line[i]));
            } else {
              mv.recordControlValue(controlCasePairings.get(sampleId), Integer.parseInt(line[i]));
            }
          } else {
            double value = Double.parseDouble(line[i]);
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
    // todo: this should be a constant from kd-match
    String expectedHeader = "id\tstatus\tmatched_case_id\n";
    BufferedReader reader = Files.getAppropriateReader(statusFile.toString());
    String actualheader = reader.readLine();

    assert actualheader.equals(expectedHeader); // todo

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
