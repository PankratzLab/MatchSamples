package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.StringJoiner;
import java.util.Vector;
import java.util.concurrent.ExecutionException;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Matrix;
import org.pankratzlab.common.PSF;
import org.pankratzlab.kdmatch.KDMatch;
import org.pankratzlab.kdmatch.KDTree;
import org.pankratzlab.kdmatch.Match;
import org.pankratzlab.kdmatch.Sample;
import org.pankratzlab.kdmatch.SelectOptimizedNeighbors;

import com.google.common.primitives.Ints;

public class MatchMaker {
  private static final String MATCH_NAIVE_TXT = "match.naive.txt";
  private static final String STATUS_NAIVE_TXT = "status.naive.txt";
  private static final String MATCH_OPTIMIZED_TXT = "match.optimized.txt";
  private static final String STATUS_OPTIMIZED_TXT = "status.optimized.txt";

  private static List<Match> kdMatchMaker(Path baseDir, Path inputSamples, List<Sample> caseList,
                                          List<Sample> controlList,
                                          HashMap<Integer, Double> numericColumnsToUseForClustering,
                                          int initialNumSelect, int finalNumSelect,
                                          FactorLoadings factorLoadings, boolean skipOptimization,
                                          int threads, Logger log) {

    KDTree<Sample> kdTree = new KDTree<>(numericColumnsToUseForClustering.keySet().size());
    log.info("Assuming 1 ID column and " + (numericColumnsToUseForClustering.keySet().size())
             + " data columns");
    log.info("Building tree for: " + caseList.get(0).getGroup());

    KDTree.addSamplesToTree(kdTree, controlList.stream());

    log.info("selecting initial " + initialNumSelect + " nearest neighbors for "
             + caseList.get(0).getGroup());

    List<Match> naiveMatches = KDTree.getNearestNeighborsForSamples(kdTree, caseList.stream(),
                                                                    initialNumSelect)
                                     .collect(Collectors.toList());

    String outputNaiveFileName = baseDir + File.separator + MATCH_NAIVE_TXT;
    String statusNaiveFileName = baseDir + File.separator + STATUS_NAIVE_TXT;
    log.info("reporting full baseline selection of " + initialNumSelect + " nearest neighbors to "
             + outputNaiveFileName);
    LinkedHashSet<String> setConvert = new LinkedHashSet<String>();

    try {
      BufferedReader headerRead = org.pankratzlab.common.Files.getAppropriateReader(inputSamples.toString());
      String[] header = headerRead.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      for (int i : numericColumnsToUseForClustering.keySet()) {
        setConvert.add(header[i]);
      }
      try {
        KDMatch.writeToFile(naiveMatches.stream(), outputNaiveFileName,
                            setConvert.stream().toArray(String[]::new),
                            setConvert.stream().toArray(String[]::new), initialNumSelect);

        KDMatch.writeSampleStatusFile(naiveMatches.stream(), statusNaiveFileName, initialNumSelect);
      } catch (IOException e) {
        e.printStackTrace();
      }
    } catch (Exception e) {
      e.printStackTrace();
      System.exit(1);
    }

    if (skipOptimization) {
      log.info("Skip optimization was flagged, returning only naive matches");
      return naiveMatches;

    } else {
      String outputOptFileName = baseDir + File.separator + MATCH_OPTIMIZED_TXT;
      String statusOptFileName = baseDir + File.separator + STATUS_OPTIMIZED_TXT;
      log.info("selecting optimized nearest neighbors");

      List<Match> optimizedMatches = null;
      try {
        optimizedMatches = SelectOptimizedNeighbors.optimizeDuplicates(naiveMatches, finalNumSelect,
                                                                       threads, log)
                                                   .collect(Collectors.toList());
        log.info("reporting optimized selection of " + finalNumSelect + " nearest neighbors to "
                 + outputOptFileName);
        KDMatch.writeToFile(optimizedMatches.stream(), outputOptFileName,
                            setConvert.stream().toArray(String[]::new),
                            setConvert.stream().toArray(String[]::new), finalNumSelect);
        KDMatch.writeSampleStatusFile(optimizedMatches.stream(), statusOptFileName, finalNumSelect);
      } catch (StackOverflowError s1) {
        s1.printStackTrace();
        log.info("To potentially prevent this StackOverflowError, try increasing the Thread Stack Size with the -Xss argument passed to the java virtual machine (i.e. java -Xss10m)");
      } catch (InterruptedException e1) {
        e1.printStackTrace();
      } catch (ExecutionException e2) {
        e2.printStackTrace();
      } catch (IOException e3) {
        e3.printStackTrace();
      }

      return optimizedMatches;
    }

  }

  private static Sample parseSample(String[] sampleLine, int idCol,
                                    HashMap<Integer, Double> numericColumnsToUseForClustering,
                                    int[] factorColumnsToAssignGroup,
                                    FactorLoadings factorLoadings) {
    StringJoiner group = new StringJoiner("_");
    String id = sampleLine[idCol];
    int status = Integer.parseInt(sampleLine[idCol + 1]);
    double[] dim = new double[numericColumnsToUseForClustering.keySet().size()];
    // for (int i = 0; i < dim.length; i++) {
    int dimIndex = 0;
    for (Entry<Integer, Double> e : numericColumnsToUseForClustering.entrySet()) {
      // TODO improve: This isn't great - I think it requires the factors input
      // argument from user
      // to be in file column order
      dim[dimIndex] = Double.parseDouble(sampleLine[e.getKey()]) * e.getValue();
      dimIndex++;
    }
    for (int i = 0; i < factorColumnsToAssignGroup.length; i++) {
      group.add(sampleLine[factorColumnsToAssignGroup[i]]);
    }
    return new Sample(id, dim, status, group.toString());
  }

  public static Path runMatching(Path dir, Path inputSamples, FactorLoadings factorLoadings,
                                 int initialNumSelect, int finalNumSelect, int threads,
                                 boolean normalize, boolean skipOptimization,
                                 Logger log) throws IOException {

    if (normalize) {
      inputSamples = normalizeFactors(dir, inputSamples, factorLoadings, log);
      log.info("Normalized input factors and wrote to file: " + inputSamples.toString());
    }

    inputSamples = handleNominalVariables(dir, inputSamples, factorLoadings);

    HashMap<Integer, Double> numericColumnsToUseForClustering = getNumericColumnsForClustering(inputSamples,
                                                                                               factorLoadings);
    int[] factorColumnsToAssignGroup = getLoadingIndices(inputSamples, factorLoadings.getFactors(),
                                                         true, true, log);
    int idColumn = 0;
    Map<String, List<Sample>> casesGroupedByStringFactor = new HashMap<String, List<Sample>>();
    Map<String, List<Sample>> controlsGroupedByStringFactor = new HashMap<String, List<Sample>>();
    Files.lines(inputSamples).map(l -> l.split("\t")).skip(1)
         .map(k -> parseSample(k, idColumn, numericColumnsToUseForClustering,
                               factorColumnsToAssignGroup, factorLoadings))
         .filter(Sample::isValidCaseOrControl).forEach(s -> {
           if (s.isCase()) {
             if (!casesGroupedByStringFactor.containsKey(s.getGroup())) {
               casesGroupedByStringFactor.put(s.getGroup(), new ArrayList<Sample>());
             }
             casesGroupedByStringFactor.get(s.getGroup()).add(s);
           } else if (s.isControl()) {
             if (!controlsGroupedByStringFactor.containsKey(s.getGroup())) {
               controlsGroupedByStringFactor.put(s.getGroup(), new ArrayList<Sample>());
             }
             controlsGroupedByStringFactor.get(s.getGroup()).add(s);
           }
         });
    final Path tempInputSamples = inputSamples;
    List<Match> results = casesGroupedByStringFactor.entrySet().stream()
                                                    .map(e -> kdMatchMaker(dir, tempInputSamples,
                                                                           e.getValue(),
                                                                           controlsGroupedByStringFactor.get(e.getKey()),
                                                                           numericColumnsToUseForClustering,
                                                                           initialNumSelect,
                                                                           finalNumSelect,
                                                                           factorLoadings,
                                                                           skipOptimization,
                                                                           threads, log))
                                                    .flatMap(List::stream)
                                                    .collect(Collectors.toList());

    return inputSamples;

  }

  /**
   * Find the index of each factor name in the sample file header and return a map from index to
   * loading (index of factor in header) -> (factor loading)
   * @param sampleFile the sample file to be read
   * @param factorloadings An object containing the factor names and loadings
   * @return A map of Integer -> Double representing index -> loading
   */
  public static HashMap<Integer, Double> getNumericColumnsForClustering(Path sampleFile,
                                                                        FactorLoadings factorloadings) {
    HashMap<Integer, Double> columnsToUse = new HashMap<>();
    try (BufferedReader origSamplesFile = org.pankratzlab.common.Files.getAppropriateReader(sampleFile.toString())) {
      Set<String> numericFactorNames = new HashSet<>(factorloadings.getNumericFactorNames());
      String[] header = origSamplesFile.readLine().strip().split(PSF.Regex.GREEDY_WHITESPACE);

      for (int j = 0; j < header.length; j++) {
        if (numericFactorNames.contains(header[j])) {
          columnsToUse.put(j, factorloadings.getLoadingForFactor(header[j]));
        }
      }

    } catch (IOException ioe) {
      ioe.printStackTrace();
      System.exit(1);
    }

    return columnsToUse;
  }

  private static Path handleNominalVariables(Path dir, Path inputSamples,
                                             FactorLoadings factorLoadings) {

    try (BufferedReader origSamplesFile = org.pankratzlab.common.Files.getAppropriateReader(inputSamples.toString())) {
      String line = origSamplesFile.readLine();

      boolean containsNomVars = false;
      ArrayList<Integer> nomInds = getNominalIndices(line, factorLoadings.getNominalFactorNames());
      if (nomInds.size() > 0) {
        containsNomVars = true;
      }
      if (!containsNomVars) {
        return inputSamples;
      }

      Path nominalized = Paths.get(dir + "/nominalized_samples.txt");

      HashMap<Integer, Set<String>> nomStruct = buildNominalStructure(inputSamples, nomInds);

      PrintWriter nomSamplesFile = org.pankratzlab.common.Files.getAppropriateWriter(nominalized.toString());

      // build header
      StringJoiner j = new StringJoiner("\t");

      String[] origLine = line.trim().split(PSF.Regex.GREEDY_WHITESPACE);
      for (int i = 0; i < origLine.length; i++) {
        if (nomStruct.keySet().contains(i)) {
          int stop = 0;
          for (String s : nomStruct.get(i)) {
            if (stop < nomStruct.get(i).size() - 1) {
              j.add(origLine[i] + "_" + s);
            }
            stop++;
          }
        } else {
          j.add(origLine[i]);
        }
      }

      nomSamplesFile.println(j.toString()); // print header to new file
      line = origSamplesFile.readLine();

      while (line != null) {
        j = new StringJoiner("\t");
        origLine = line.trim().split(PSF.Regex.GREEDY_WHITESPACE);
        for (int i = 0; i < origLine.length; i++) {
          if (nomStruct.keySet().contains(i)) {
            int stop = 0;
            for (String s : nomStruct.get(i)) {
              if (stop < nomStruct.get(i).size() - 1) {
                if (s.equalsIgnoreCase(origLine[i])) {
                  j.add("1");
                } else {
                  j.add("0");
                }
              }
              stop++;
            }
          } else {
            j.add(origLine[i]);
          }
        }
        nomSamplesFile.println(j.toString());
        line = origSamplesFile.readLine();
      }
      nomSamplesFile.close();
      return nominalized;
    } catch (IOException e) {
      System.out.println("Sample file " + inputSamples + " not found.");
      System.exit(1);
    }
    return null;

  }

  private static ArrayList<Integer> getNominalIndices(String header,
                                                      ArrayList<String> nominalFactorNames) {
    String[] split = header.trim().split(PSF.Regex.GREEDY_WHITESPACE);
    ArrayList<Integer> indices = new ArrayList<Integer>();
    for (int i = 0; i < split.length; i++) {
      if (nominalFactorNames.contains(split[i])) {
        indices.add(i);
      }
    }
    return indices;
  }

  private static HashMap<Integer, Set<String>> buildNominalStructure(Path sampleFile,
                                                                     ArrayList<Integer> nominalIndices) {
    HashMap<Integer, Set<String>> noms = new HashMap<Integer, Set<String>>();
    for (Integer x : nominalIndices) {
      noms.put(x, new LinkedHashSet<String>());
    }
    try (BufferedReader uniqueNomFinder = org.pankratzlab.common.Files.getAppropriateReader(sampleFile.toString())) {
      String header = uniqueNomFinder.readLine(); // skip header
      String lineString = uniqueNomFinder.readLine();
      while (lineString != null) {
        String[] line = lineString.trim().split(PSF.Regex.GREEDY_WHITESPACE);
        for (Integer x : nominalIndices) {
          noms.get(x).add(line[x]);
        }
        lineString = uniqueNomFinder.readLine();
      }
    } catch (IOException e) {
      System.out.println("Sample file " + sampleFile + " not found.");
      System.exit(1);
    }
    return noms;
  }

  private static Path normalizeFactors(Path dir, Path inputSamples, FactorLoadings factors,
                                       Logger log) throws IOException {
    String normalizedSamples = dir + "/normalized.txt";
    String[][] matrix;
    double[][] allData;
    int[] factorIndices = getLoadingIndices(inputSamples, factors.getFactors(), false, false, log);
    matrix = HashVec.loadFileToStringMatrix(inputSamples.toString(), true, factorIndices,
                                            PSF.Regex.GREEDY_WHITESPACE, 1000, false);

    allData = new double[factorIndices.length][];
    for (int i = 0; i < factorIndices.length; i++) {
      allData[i] = ArrayUtils.toDoubleArray(Matrix.extractColumn(matrix, i));
      allData[i] = ArrayUtils.normalize(allData[i]);
    }
    BufferedReader inputReader = org.pankratzlab.common.Files.getAppropriateReader(inputSamples.toString());
    PrintWriter normResults = org.pankratzlab.common.Files.getAppropriateWriter(normalizedSamples);
    String inputLine = inputReader.readLine();
    String[] header = inputLine.trim().split(PSF.Regex.GREEDY_WHITESPACE);
    normResults.println(inputLine); // header
    String[] line;
    int row = 0;
    inputLine = inputReader.readLine(); // go past header
    while (inputLine != null) {
      line = inputLine.trim().split(PSF.Regex.GREEDY_WHITESPACE);
      StringJoiner join = new StringJoiner("\t");
      Boolean normalizedColumn = false;
      int normColIndex = 0;
      for (int i = 0; i < header.length; i++) {
        for (int f : factorIndices) {
          if (f == i) {
            normalizedColumn = true;
            normColIndex = Ints.indexOf(factorIndices, i);
          }
        }
        if (normalizedColumn) {
          join.add(allData[normColIndex][row] + "");
        } else {
          join.add(line[i]);
        }
        normalizedColumn = false;
      }
      row++;
      normResults.println(join);
      inputLine = inputReader.readLine();
    }
    normResults.close();
    return Paths.get(normalizedSamples);

  }

  private static int[] getLoadingIndices(Path sampleFile, Map<String, String> factorLoadings,
                                         boolean force, boolean includeNominal, Logger log) {
    Set<String> forcedCols = new LinkedHashSet<String>();
    Set<String> doubleCols = new LinkedHashSet<String>();
    boolean isNom;
    double temp;
    for (Entry<String, String> s : factorLoadings.entrySet()) {
      isNom = false;
      if (s.getValue().equalsIgnoreCase("force")) {
        forcedCols.add(s.getKey());
      } else {
        try {
          if (s.getValue().contains("nom_")) {
            temp = Double.parseDouble(s.getValue().trim().split("_")[1]);
            isNom = true;
          } else {
            temp = Double.parseDouble(s.getValue());
          }
          if (includeNominal) {
            doubleCols.add(s.getKey());
          } else if (!includeNominal && !isNom) {
            doubleCols.add(s.getKey());
          } else {
            continue;
          }
        } catch (NumberFormatException nfe) {
          log.info("Invalid loading found: " + s);
          System.exit(1);
        }
      }
    }

    int[] loadingIndices;
    if (force) {
      loadingIndices = new int[forcedCols.size()];
    } else {
      loadingIndices = new int[doubleCols.size()];
    }

    int currentIndex = 0;
    try (BufferedReader r = org.pankratzlab.common.Files.getAppropriateReader(sampleFile.toAbsolutePath()
                                                                                        .toString())) {
      String[] header = r.readLine().split("\t");
      for (int i = 0; i < header.length; i++) {
        if (force) {
          if (forcedCols.contains(header[i])) {
            loadingIndices[currentIndex] = i;
            currentIndex++;
          }
        } else {
          if (doubleCols.contains(header[i])) {
            loadingIndices[currentIndex] = i;
            currentIndex++;
          }
        }
      }
    } catch (IOException ioe) {
      log.info("Can't find samples file: " + sampleFile);
    }
    return loadingIndices;
  }

  public static String buildVisHelpers(Path dir, Path fullResultsFile, Path samplesFile,
                                       int currentIteration, Logger log) throws IOException {
    Path visDir = Paths.get(dir + "/visual_helpers/");
    if (!visDir.toFile().exists()) {
      Files.createDirectory(visDir);
    }

    Path tempVisHelperFile = Paths.get(visDir + "/vis_helper_" + (currentIteration + 1) + ".temp");
    Path tempFactorsHelperFile = Paths.get(visDir + "/vis_helper_factors.temp");

    File f = tempVisHelperFile.toFile();
    File f2 = tempFactorsHelperFile.toFile();
    if (f.exists()) {
      f.delete();
    }
    if (f2.exists()) {
      f2.delete();
    }

    PrintWriter factorsHelperWriter;
    Vector<String> factorsData = HashVec.loadFileToVec(samplesFile.toString(), false, false, true);
    factorsHelperWriter = org.pankratzlab.common.Files.getAppropriateWriter(tempFactorsHelperFile.toString());
    String[] factorsLine = factorsData.elementAt(0).split(PSF.Regex.GREEDY_WHITESPACE);
    StringJoiner header = new StringJoiner("\t");
    for (int i = 0; i < factorsLine.length; i++) {
      if (i != 1) {
        header.add(factorsLine[i]);
      }
    }
    factorsHelperWriter.println(header);
    for (int i = 1; i < factorsData.size(); i++) {
      factorsLine = factorsData.elementAt(i).trim().split(PSF.Regex.GREEDY_WHITESPACE);
      StringJoiner newLine = new StringJoiner("\t");
      for (int s = 0; s < factorsLine.length; s++) {
        if (s != 1) {
          newLine.add(factorsLine[s]);
        }
      }
      factorsHelperWriter.println(newLine);
    }
    factorsHelperWriter.close();

    PrintWriter visHelperWriter;
    int controlColumn = -1;
    Vector<String> data = HashVec.loadFileToVec(fullResultsFile.toString(), false, false, true);
    visHelperWriter = org.pankratzlab.common.Files.getAppropriateWriter(tempVisHelperFile.toString());
    String[] line;
    visHelperWriter.println("Case\tControl\tDistance");
    line = data.elementAt(0).trim().split(PSF.Regex.GREEDY_WHITESPACE);
    int index = 0;
    for (String s : line) {
      if (s.equalsIgnoreCase("control_" + (currentIteration + 1) + "_id")) {
        controlColumn = index;
        break;
      }
      index++;
    }
    for (int i = 1; i < data.size(); i++) {
      line = data.elementAt(i).trim().split(PSF.Regex.GREEDY_WHITESPACE);
      visHelperWriter.println(line[0] + "\t" + line[controlColumn] + "\t"
                              + line[controlColumn + 1]);
    }
    visHelperWriter.close();
    return tempVisHelperFile.toString();
  }

  public static List<String> parseEvalArgs(String argString) {
    List<String> evalArgs = new ArrayList<>();
    String[] commaSplit = argString.split(",");

    for (String s : commaSplit) {
      evalArgs.add(s.split(":")[0]);
    }
    return evalArgs;
  }

  public static void main(String[] args) {

    Path d = Paths.get("./");
    FactorLoadings factorLoadings = null;
    Path samples = Paths.get(d + "samples.txt");
    boolean normalize = true;
    int numArgs = args.length;
    int finalNumSelect = 4;
    int multiplier = 5;
    int threads = Runtime.getRuntime().availableProcessors();
    boolean vis = false;
    boolean onlyBuildVisFiles = false;
    boolean skipEval = false;
    boolean skipOptimization = false;
    List<String> evalArgs = null;
    Logger log;

    String usage = "\n" + "gwas.MatchMaker requires at least 1 argument\n"
                   + "(1) Working directory (e.g. dir=./ (default))\n"
                   + "(2) Samples filename with ID column, case/control column (1/0), and factor columns (e.g. samples=samples.txt (default))\n"
                   + "(3) Factors columns to use and their loadings; *keep in file order (e.g. factors=PC1:4,PC2:4,sex:force,phenograph:force)\n"
                   + "(4) Normalize inputs (e.g. normalize=true (default))\n"
                   + "(5) Iterations - number of controls to match to each case (e.g. iterations=4 (default))\n"
                   + "(6) Multiplier - naive matching multiplier (e.g. multiplier=5 (default))\n"
                   + "(7) Normalize the input factors before matching (e.g. normalize=true (default))\n"
                   + "(8) Eval - which arguments you want to check for concordance (e.g. eval=age,sex,PC1 (defaults to all factor names))\n"
                   + "(9) Visualize results - (e.g. vis=false (default))\n"
                   + "(10) Only build the visualizer files to run separately - (e.g. onlyBuildVisFiles=false (default))\n"
                   + "(11) Skip evaluation - skip generating statistical analysis of matchmaking quality (e.g. skipEval=false (default))\n"
                   + "(12) Skip optimization - skip optimization of matches (e.g. skipOptimization=false (default))\n"
                   + "(13) Number of threads to use (default = Runtime.getRuntime().availableProcessors()) (e.g. threads=10) ";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        d = Paths.get(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("samples=")) {
        samples = Paths.get(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("factors=")) {
        factorLoadings = new FactorLoadings(splitEq(arg));
        numArgs--;
      } else if (arg.startsWith("normalize=")) {
        normalize = Boolean.parseBoolean(splitEq(arg));
      } else if (arg.startsWith("iterations=")) {
        finalNumSelect = Integer.parseInt(splitEq(arg));
      } else if (arg.startsWith("multiplier=")) {
        multiplier = Integer.parseInt(splitEq(arg));
      } else if (arg.startsWith("eval=")) {
        evalArgs = parseEvalArgs(splitEq(arg));
      } else if (arg.startsWith("vis=")) {
        vis = Boolean.parseBoolean(splitEq(arg));
      } else if (arg.startsWith("onlyBuildVisFiles=")) {
        onlyBuildVisFiles = Boolean.parseBoolean(splitEq(arg));
      } else if (arg.startsWith("skipEval=")) {
        skipEval = Boolean.parseBoolean(splitEq(arg));
      } else if (arg.startsWith("skipOptimization=")) {
        skipOptimization = Boolean.parseBoolean(splitEq(arg));
      } else if (arg.startsWith("threads=")) {
        threads = Integer.parseInt(splitEq(arg));
      }
    }

    if (factorLoadings == null) {
      System.err.println("No factors supplied, unable to continue");
      System.err.println(usage);
      System.exit(1);
    }

    int initialNumSelect = finalNumSelect * multiplier;
    samples = Paths.get(d + File.separator + samples);

    log = Logger.getAnonymousLogger();
    log.info("Starting sample match using k-d tree nearest neighbors.");

    try {
      List<String> fileNames = List.of(MATCH_NAIVE_TXT, STATUS_NAIVE_TXT, MATCH_OPTIMIZED_TXT,
                                       STATUS_OPTIMIZED_TXT);
      final Path finalD = d;
      boolean outputExists = fileNames.stream()
                                      .map(name -> new File(finalD + File.separator + name))
                                      .anyMatch(File::exists);
      if (outputExists) {
        log.info("Output already exists.");
        System.exit(0);
      }
      Path normalizedSamples = runMatching(d, samples, factorLoadings, initialNumSelect,
                                           finalNumSelect, threads, normalize, skipOptimization,
                                           log);
      if (vis) {
        HashMap<Integer, Double> temp = getNumericColumnsForClustering(normalizedSamples,
                                                                       factorLoadings);
        int[] loadingIndicesForVis = new int[temp.keySet().size()];
        int ind = 0;
        for (Integer x : temp.keySet()) {
          loadingIndicesForVis[ind] = x;
          ind++;
        }
        for (int s = 0; s < loadingIndicesForVis.length; s++) {
          loadingIndicesForVis[s] = loadingIndicesForVis[s] - 1;
        }
        for (int i = 0; i < finalNumSelect; i++) {
          Path resultsFile = Paths.get(new File(d + File.separator
                                                + MATCH_OPTIMIZED_TXT).toString());
          buildVisHelpers(d, resultsFile, normalizedSamples, i, log);
          if (!onlyBuildVisFiles) {
            new MatchesVisualized(d.toString(), normalizedSamples.toString(),
                                  d + "/visual_helpers/vis_helper_factors.temp",
                                  loadingIndicesForVis,
                                  d + "/visual_helpers/vis_helper_" + (i + 1) + ".temp", true);
          }
        }
        // FileUtils.deleteDirectory(new File(d + "/visual_helpers/"));
      }

      if (!skipEval) {
        log.info("Performing eval");
        MatchingVariable[] matchingVariables;
        if (evalArgs == null) {
          matchingVariables = MatchingVariable.fromNames(factorLoadings.getFactors().keySet());
        } else {
          matchingVariables = MatchingVariable.fromNames(evalArgs);
        }

        File naiveStatusFile = new File(d + File.separator + STATUS_NAIVE_TXT);
        REval rEvalNaive = new REval(matchingVariables, naiveStatusFile, samples.toFile());
        rEvalNaive.readPhenotypeFile();
        File naiveRevalOutputFile = new File(d + File.separator + "eval_results_naive.tsv");
        rEvalNaive.writeTableOutputToFile(naiveRevalOutputFile);
        log.info("Wrote eval output to " + naiveRevalOutputFile.getName());

        if (skipOptimization) {
          log.info("Skip optimization was flagged, skipping evaluation of optimized matches");
        } else {
          File optimizedStatusFile = new File(d + File.separator + STATUS_OPTIMIZED_TXT);
          REval rEvalOptimized = new REval(matchingVariables, optimizedStatusFile,
                                           samples.toFile());
          rEvalOptimized.readPhenotypeFile();
          File optimizedRevalOutputFile = new File(d + File.separator
                                                   + "eval_results_optimized.tsv");
          rEvalOptimized.writeTableOutputToFile(optimizedRevalOutputFile);
          log.info("Wrote eval output to " + optimizedRevalOutputFile.getName());
        }
        // obviously there's an opportunity for optimization here because this is reading the entire
        // phenotype file twice, when only the statuses have changed.
      }

    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  private static String splitEq(String arg) {
    return arg.split("=")[1];
  }

}
