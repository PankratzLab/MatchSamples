package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.pankratzlab.common.Files;
import org.pankratzlab.common.PSF;

public class REval {
  final MatchingVariable[] matchingVariables;
  final File statusFile;
  final File phenotypeFile;
  final Map<String, String> pairings;

  public REval(MatchingVariable[] matchingVariables, File statusFile,
               File phenotypeFile) throws IOException {
    this.matchingVariables = matchingVariables;
    // todo: don't use assert
    assert statusFile.exists();
    assert phenotypeFile.exists();

    this.statusFile = statusFile;
    this.phenotypeFile = phenotypeFile;

    this.pairings = readPairings();
  }

  public void evaluate() throws IOException {
    BufferedReader phenoReader = Files.getAppropriateReader(phenotypeFile.toString());
    String[] phenoHeader = phenoReader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);

    Map<Integer, MatchingVariable> matchingVarIndexMap = new HashMap<>();
    for (MatchingVariable mv : matchingVariables) {
      int i = mv.findIndexInHeader(phenoHeader);
      assert i != -1; // todo: check that all the matching variables end up in the map
      matchingVarIndexMap.put(i, mv);
      // todo: there should be some check that no two matching variables share an index
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
      line = inputLine.trim().split(PSF.Regex.GREEDY_WHITESPACE);
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
}
