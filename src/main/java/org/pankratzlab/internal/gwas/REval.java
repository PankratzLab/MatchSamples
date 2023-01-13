package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.pankratzlab.common.Files;
import org.pankratzlab.common.PSF;
import org.pankratzlab.kdmatch.KDMatch;

import static java.lang.System.exit;

public class REval {
  final MatchingVariable[] matchingVariables;
  final File statusFile;
  final File phenotypeFile;
  final Map<String, String> controlCasePairings;
  private boolean haveReadPhenotype = false;

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

  public void readPhenotypeFile() throws IOException {
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
    haveReadPhenotype = true;
  }

  public void writeTableOutputToFile(File outputFile) throws IOException {
    if (!haveReadPhenotype) {
      readPhenotypeFile();
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

  public static void main(String[] args) throws IOException {

    //@format:off
    String usage = "\ngwas.REval usage:\n"
                   + "dir=working_directory/ (optional, default=./)\n"
                   + "status=path/to/status_file.tsv (optional, default=status.optimized.txt, path is relative to working directory)\n"
                   + "phenotype=path/to/phenotype_file.tsv (required, path is relative to working directory)\n"
                   + "matchingVars=foo;bar (required, a semicolon-separated list of column names in the phenotype file)"
                   + "output=path/to/output_file.tsv (optional, default=reval_results.tsv, path is relative to working directory\n"
                   + "\n"
                   + "status file:        This should be the output from MatchMaker\n"
                   + "phenotype file:     This should be the phenotype file that was used as input to produce the status file\n"
                   + "matching vars file: Two column tsv describing the matching variables to be evaluated\n"
                   + "    header_name -> the name of a column in the phenotype file\n"
                   + "    is_binary   -> true or false, if true, treat this variable as having only two possible values"
                   + "\n";
    //@format:on
    Path dir = Paths.get("./");
    Path status = Paths.get("status.optimized.txt");
    Path phenotype = null;
    String matchingVariables = null;
    Path output = Paths.get("reval_results.tsv");

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = parsePathFromArg(arg);
      } else if (arg.startsWith("status=")) {
        status = parsePathFromArg(arg);
      } else if (arg.startsWith("phenotype=")) {
        phenotype = parsePathFromArg(arg);
      } else if (arg.startsWith("matchingVars=")) {
        matchingVariables = arg.split("=")[1];
      } else if (arg.startsWith("output=")) {
        output = parsePathFromArg(arg);
      }
    }

    if (phenotype == null || matchingVariables == null) {
      System.err.println("No phenotype file provided, cannot continue.");
      System.err.println(usage);
      exit(1);
    }

    MatchingVariable[] matchingVariablesArray = MatchingVariable.fromSemicolonSeparatedString(matchingVariables);
    File statusFile = dir.resolve(status).toFile();
    File phenotypeFile = dir.resolve(phenotype).toFile();
    File outputFile = dir.resolve(output).toFile();

    if (!outputFile.getParentFile().isDirectory()) {
      throw new IllegalArgumentException("The directory containing the specified output file does not exist.");
    }
    if (outputFile.isFile()) {
      throw new IllegalArgumentException("The specified output file already exists.");
    }

    REval rEval = new REval(matchingVariablesArray, statusFile, phenotypeFile);
    rEval.readPhenotypeFile();
    rEval.writeTableOutputToFile(outputFile);
  }

  private static Path parsePathFromArg(String arg) {
    return Paths.get(arg.split("=")[1]);
  }
}
