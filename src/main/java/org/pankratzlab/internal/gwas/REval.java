package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
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
  private boolean evaluated = false;

  private final DataBox dataBox;

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

    if (!statusFile.isFile()) {
      throw new IllegalArgumentException("Provided status file does not exist or is not a normal file.");
    }
    if (!phenotypeFile.isFile()) {
      throw new IllegalArgumentException("Provided phenotype file does not exist or is not a normal file.");
    }

    this.statusFile = statusFile;
    this.phenotypeFile = phenotypeFile;

    this.controlCasePairings = readPairings();

    this.dataBox = new DataBox(matchingVariables, controlCasePairings);

  }

  public void evaluate() throws IOException {
    BufferedReader phenoReader = Files.getAppropriateReader(phenotypeFile.toString());
    String[] phenoHeader = splitTsvLine(phenoReader.readLine());

    // figure out which column corresponds to each matching variable
    for (MatchingVariable mv : matchingVariables) {
      int i = mv.findIndexInHeader(phenoHeader);
      if (i == -1) {
        throw new IllegalStateException("Unable to find column for matching variable; expected header name: "
                                        + mv.headerName);
      }
    }

    String[] line;
    String inputLine = phenoReader.readLine();
    while (inputLine != null) {
      line = splitTsvLine(inputLine);
      dataBox.recordData(line);
      inputLine = phenoReader.readLine();
    }
    evaluated = true;
  }

  public void writeTableOutputToFile(File outputFile) throws IOException {
    if (!evaluated) {
      evaluate();
    }
    if (outputFile.exists()) {
      throw new IllegalArgumentException("Provided output file already exists!");
    }
    String header = "variable_name\tcase_avg\tcontrol_avg\tconcordance\tunivariate_p\tmultivariate_p";

    try (PrintWriter writer = new PrintWriter(outputFile)) {
      writer.println(header);
      for (MatchingVariable mv : matchingVariables) {
        writer.println(mv.getTableLine());
      }
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
      inputLine = reader.readLine();
    }
    return pairings;
  }

  private static String[] splitTsvLine(String inputLine) {
    return inputLine.trim().split(PSF.Regex.GREEDY_WHITESPACE);
  }
}
