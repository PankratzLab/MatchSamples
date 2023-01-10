package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
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
  private int totalSampleCount = 0;

  private final BinaryDataBox binaryDataBox;
  private final ContinuousDataBox continuousDataBox;

  public REval(MatchingVariable[] matchingVariables, File statusFile,
               File phenotypeFile) throws IOException {
    Set<String> headerNamesSeen = new HashSet<>();

    for (MatchingVariable mv : matchingVariables) {
      if (headerNamesSeen.contains(mv.headerName)) {
        throw new IllegalArgumentException("Multiple matching variables provided for this header name: "
                                           + mv.headerName);
      } else {
        headerNamesSeen.add(mv.headerName);
      }
    }
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

    MatchingVariable[] binaryMvs = Arrays.stream(matchingVariables).filter(mv -> mv.isBinary)
                                         .toArray(MatchingVariable[]::new);
    binaryDataBox = new BinaryDataBox(binaryMvs, totalSampleCount, controlCasePairings, cases);

    MatchingVariable[] continuousMvs = Arrays.stream(matchingVariables).filter(mv -> !mv.isBinary)
                                             .toArray(MatchingVariable[]::new);
    continuousDataBox = new ContinuousDataBox(continuousMvs, totalSampleCount, controlCasePairings,
                                              cases);

  }

  public void evaluate() throws IOException {
    BufferedReader phenoReader = Files.getAppropriateReader(phenotypeFile.toString());
    String[] phenoHeader = splitTsvLine(phenoReader.readLine());

    // figure out which column corresponds to each matching variable
    Map<Integer, MatchingVariable> matchingVarsByIndex = new HashMap<>();
    for (MatchingVariable mv : matchingVariables) {
      int i = mv.findIndexInHeader(phenoHeader);
      if (i == -1) {
        throw new IllegalStateException("Unable to find column for matching variable; expected header name: "
                                        + mv.headerName);
      }
      matchingVarsByIndex.put(i, mv);
    }

    String[] line;
    String inputLine = phenoReader.readLine();
    while (inputLine != null) {
      line = splitTsvLine(inputLine);
      String sampleId = line[0];
      if (cases.contains(sampleId) || controls.contains(sampleId)) {
        for (int matchVarIndex : matchingVarsByIndex.keySet()) {
          MatchingVariable mv = matchingVarsByIndex.get(matchVarIndex);
          if (mv.isBinary) {
            binaryDataBox.recordValue(sampleId, mv, Integer.parseInt(line[matchVarIndex]));
          } else {
            double value = Double.parseDouble(line[matchVarIndex]);
            continuousDataBox.recordValue(mv, sampleId, value);
          }
        }
      }
      inputLine = phenoReader.readLine();
    }
  }

  private void univariateLeastSquares() {

    // LeastSquares leastSquares = new LeastSquares.LeastSquaresBuilder();
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
      totalSampleCount++;
      inputLine = reader.readLine();
    }
    return pairings;
  }

  private static String[] splitTsvLine(String inputLine) {
    return inputLine.trim().split(PSF.Regex.GREEDY_WHITESPACE);
  }
}
