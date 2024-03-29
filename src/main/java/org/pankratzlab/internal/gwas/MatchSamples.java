package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import org.pankratzlab.kdmatch.*;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.DoubleVector;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.Matrix;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.mining.Distance;
import org.pankratzlab.common.mining.Transformations;
import org.pankratzlab.common.stats.Correlation;
import org.pankratzlab.common.stats.Ttest;
import org.pankratzlab.utils.gwas.Plink;

import com.google.common.primitives.Doubles;

public class MatchSamples {

  private String dir;
  private String samplesFile;
  private String factorsFile;
  private String minOrMax;
  private String[] cols;
  private boolean skipAge;
  private double[] loadings;
  private Logger log;

  public MatchSamples(String directory, String samplesFile, String factors, String[] cols,
                      double[] loadings, String minOrMax, boolean skipVisualizer,
                      boolean skipAgeSexMatch) {
    this.dir = directory;
    this.samplesFile = dir + samplesFile;
    this.factorsFile = dir + factors;
    this.loadings = loadings;
    this.cols = cols;
    this.minOrMax = minOrMax;
    this.skipAge = skipAgeSexMatch;
  }

  /*
   * public static String kdMatchMaker(Path dir, Path cases, Path controls, String[] factorTargets,
   * double[] factorLoadings, boolean normalizeFactors, Logger log) throws IOException {
   * 
   * String[] header = java.nio.file.Files.lines(factors).findFirst().get().toString().trim().split(
   * "\t"); KDTree<Sample> kdTree = new KDTree<>(header.length - 1);
   * 
   * log.report("Assuming 1 ID column and " + (header.length - 1) + " data columns");
   * log.report("building tree from " + factors.toString());
   * 
   * KDTree.addSamplesToTree(kdTree, KDMatch.getSampleStreamFromFile(factors));
   * 
   * log.report("selecting initial " + initialNumSelect + " nearest neighbors for " +
   * inputFileAnchor.toString()); }
   */

  public static String matchMaker(String dir, String samplesFile, String factorfile,
                                  String[] factorTargets, double[] factorLoadings,
                                  boolean normalizeFactors) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, cases, controls;
    double[][] dists, newDists = null;
    double[][] allData, caseData, controlData;
    long time;
    int[] factorIndices = null;
    String filename;
    int iAnch, iBarn;
    String[] ids;
    String[][] matrix;

    if (factorTargets.length != factorLoadings.length) {
      System.err.println("Error - mismatch in the number of factorTargets/factorLoadings");
      System.exit(1);
    }

    filename = "distances_";
    for (int i = 0; i < factorTargets.length; i++) {
      filename += (i == 0 ? "" : ",") + factorTargets[i] + "x"
                  + ext.formDeci(factorLoadings[i], 10);
    }
    filename += ".xln";

    if (!new File(dir + filename).exists()) {
      System.out.println("Creating " + filename);
      time = new Date().getTime();
      cases = samplesFileToStringArray(dir + samplesFile, dir + factorfile, 1);
      caseData = new double[cases.length][];
      controls = samplesFileToStringArray(dir + samplesFile, dir + factorfile, 0);
      controlData = new double[controls.length][];

      factorIndices = ext.indexFactors(factorTargets,
                                       Files.getHeaderOfFile(dir + factorfile,
                                                             PSF.Regex.GREEDY_WHITESPACE,
                                                             new Logger()),
                                       false);
      ids = HashVec.loadFileToStringArray(dir + factorfile, true, new int[] {0}, false);
      matrix = HashVec.loadFileToStringMatrix(dir + factorfile, true, factorIndices,
                                              PSF.Regex.GREEDY_WHITESPACE, 1000, false);
      allData = new double[factorIndices.length][];
      for (int i = 0; i < factorTargets.length; i++) {
        allData[i] = ArrayUtils.toDoubleArray(Matrix.extractColumn(matrix, i));
        if (normalizeFactors) {
          allData[i] = ArrayUtils.normalize(allData[i]);
        }
      }
      for (int i = 0; i < ids.length; i++) {
        iAnch = ext.indexOfStr(ids[i], cases);
        iBarn = ext.indexOfStr(ids[i], controls);
        if (iAnch >= 0) {
          caseData[iAnch] = new double[factorIndices.length];
          for (int j = 0; j < factorIndices.length; j++) {
            caseData[iAnch][j] = allData[j][i] * factorLoadings[j];
          }
        } else if (iBarn >= 0) {
          controlData[iBarn] = new double[factorIndices.length];
          for (int j = 0; j < factorIndices.length; j++) {
            controlData[iBarn][j] = allData[j][i] * factorLoadings[j];
          }
        }
      }

      for (int i = 0; i < cases.length; i++) {
        if (caseData[i] == null) {
          System.err.println("Error - data for anchor '" + cases[i] + "' not found in "
                             + factorfile);
        }
      }
      for (int i = 0; i < controls.length; i++) {
        if (controlData[i] == null) {
          System.err.println("Error - data for barnacle '" + controls[i] + "' not found in "
                             + factorfile);
        }
      }

      System.out.println("Initialized in " + ext.getTimeElapsed(time));
      time = new Date().getTime();
      dists = Matrix.doubleMatrix(cases.length, controls.length, -999);
      for (int i = 0; i < cases.length; i++) {
        for (int j = 0; j < controls.length; j++) {
          dists[i][j] = Distance.euclidean(caseData[i], controlData[j]);
        }
      }
      System.out.println("Finished euclidean calculations in " + ext.getTimeElapsed(time));
      time = new Date().getTime();

      try {
        writer = Files.openAppropriateWriter(dir + filename);
        // writer = new PrintWriter(new
        // FileWriter(dir+"distances_1-100.xln"));
        writer.println(cases.length + "\t" + controls.length);
        writer.println("Anchor\t" + ArrayUtils.toStr(controls));
        for (int i = 0; i < cases.length; i++) {
          writer.println(cases[i] + "\t" + ArrayUtils.toStr(dists[i]));
        }
        writer.close();
      } catch (Exception e) {
        System.err.println("Error writing distances");
        e.printStackTrace();
      }
      System.out.println("Finished writing distances_" + ArrayUtils.toStr(factorIndices, ",")
                         + " in " + ext.getTimeElapsed(time));
    } else {
      time = new Date().getTime();

      try {
        reader = new BufferedReader(new FileReader(dir + filename));
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        cases = new String[Integer.parseInt(line[0])];
        controls = new String[Integer.parseInt(line[1])];
        newDists = Matrix.doubleMatrix(cases.length, controls.length, -999);
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        for (int i = 0; i < controls.length; i++) {
          controls[i] = line[i + 1];
        }
        for (int i = 0; i < cases.length; i++) {
          line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
          cases[i] = line[0];
          newDists[i] = ArrayUtils.toDoubleArray(ArrayUtils.subArray(line, 1));
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + dir + filename + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + dir + filename + "\"");
        System.exit(2);
      }
      System.out.println("Finished reading in " + ext.rootOf(filename) + " in "
                         + ext.getTimeElapsed(time));
    }
    time = new Date().getTime();

    return filename;
  }

  public static String[] samplesFileToStringArray(String samplesFile, String factorFile,
                                                  int caseOrControl) {
    Set<String> validSamples = new HashSet<String>();
    try (BufferedReader r = Files.getAppropriateReader(factorFile)) {
      String id = r.readLine();
      while (id != null) {
        validSamples.add(id.trim().split(PSF.Regex.GREEDY_WHITESPACE)[0]);
        id = r.readLine();
      }
    } catch (IOException e) {
      System.out.println("Couldn't find factors file: " + factorFile);
      e.printStackTrace();
    }

    Vector<String> v = new Vector<String>();
    try (BufferedReader reader = Files.getAppropriateReader(samplesFile)) {
      String line = reader.readLine();
      String sample;
      if (line.contains("status")) {
        line = reader.readLine();
      }
      while (line != null) {
        sample = line.trim().split(PSF.Regex.GREEDY_WHITESPACE)[0];
        if ((Integer.parseInt(line.trim().split(PSF.Regex.GREEDY_WHITESPACE)[1]) == caseOrControl)
            && validSamples.contains(sample)) {
          v.add(line.trim().split(PSF.Regex.GREEDY_WHITESPACE)[0]);
        }
        line = reader.readLine();
      }
    } catch (IOException ioe) {
      System.out.println("Couldn't find samples file: " + samplesFile);
      ioe.printStackTrace();
    }
    return v.toArray(new String[0]);
  }

  public static void parseClusterfile(String dir, String anchorList, String barnacleList,
                                      String clusterfile) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, anchors, barnacles;
    int iAnch, iBarn;
    long time;
    double[][] pihats, dsts, ratios;

    time = new Date().getTime();
    anchors = HashVec.loadFileToStringArray(dir + anchorList, false, new int[] {0}, true);
    barnacles = HashVec.loadFileToStringArray(dir + barnacleList, false, new int[] {0}, true);
    pihats = Matrix.doubleMatrix(anchors.length, barnacles.length, -999);
    dsts = Matrix.doubleMatrix(anchors.length, barnacles.length, -999);
    ratios = Matrix.doubleMatrix(anchors.length, barnacles.length, -999);
    try {
      reader = new BufferedReader(new FileReader(dir + clusterfile));
      ext.checkHeader(reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE),
                      Plink.CLUSTER_HEADER, true);
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        iAnch = Math.max(ext.indexOfStr(line[1], anchors), ext.indexOfStr(line[3], anchors));
        iBarn = Math.max(ext.indexOfStr(line[1], barnacles), ext.indexOfStr(line[3], barnacles));
        if (iAnch != -1 && iBarn != -1) {
          pihats[iAnch][iBarn] = 1 - Double.parseDouble(line[7]);
          dsts[iAnch][iBarn] = 1 - Double.parseDouble(line[12]);
          ratios[iAnch][iBarn] = 3 - Double.parseDouble(line[16]);
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + clusterfile
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + clusterfile + "\"");
      System.exit(2);
    }
    System.out.println("Finished parsing file in " + ext.getTimeElapsed(time));
    time = new Date().getTime();

    try {
      writer = Files.openAppropriateWriter(dir + "pihats.xln");
      writer.println(anchors.length + "\t" + barnacles.length);
      writer.println("Anchor\t" + ArrayUtils.toStr(barnacles));
      for (int i = 0; i < anchors.length; i++) {
        writer.println(anchors[i] + "\t" + ArrayUtils.toStr(pihats[i]));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing pihats");
      e.printStackTrace();
    }
    try {
      writer = Files.openAppropriateWriter(dir + "dsts.xln");
      writer.println(anchors.length + "\t" + barnacles.length);
      writer.println("Anchor\t" + ArrayUtils.toStr(barnacles));
      for (int i = 0; i < anchors.length; i++) {
        writer.println(anchors[i] + "\t" + ArrayUtils.toStr(dsts[i]));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing dsts");
      e.printStackTrace();
    }
    try {
      writer = Files.openAppropriateWriter(dir + "ratios.xln");
      writer.println(anchors.length + "\t" + barnacles.length);
      writer.println("Anchor\t" + ArrayUtils.toStr(barnacles));
      for (int i = 0; i < anchors.length; i++) {
        writer.println(anchors[i] + "\t" + ArrayUtils.toStr(ratios[i]));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing ratios");
      e.printStackTrace();
    }
    System.out.println("Finished writing values in " + ext.getTimeElapsed(time));

  }

  public static void correlate(String dir, String file1, String file2) {
    BufferedReader reader;
    String[] line, anchors, barnacles;
    double[][] values = new double[2][];

    try {
      reader = new BufferedReader(new FileReader(dir + file1));
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      anchors = new String[Integer.parseInt(line[0])];
      barnacles = new String[Integer.parseInt(line[1])];
      values[0] = ArrayUtils.doubleArray(anchors.length * barnacles.length, -999);
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      for (int i = 0; i < barnacles.length; i++) {
        barnacles[i] = line[i + 1];
      }
      for (int i = 0; i < anchors.length; i++) {
        anchors[i] = line[0];
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        for (int j = 0; j < barnacles.length; j++) {
          values[0][i * barnacles.length + j] = Double.parseDouble(line[j + 1]);
        }
      }
      reader.close();

      reader = new BufferedReader(new FileReader(dir + file2));
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      if (Integer.parseInt(line[0]) != anchors.length
          || Integer.parseInt(line[1]) != barnacles.length) {
        System.err.println("Error - the two files have different numbers of anchors and barnacles");
        System.exit(1);
      }
      values[1] = ArrayUtils.doubleArray(anchors.length * barnacles.length, -999);
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      for (int i = 0; i < barnacles.length; i++) {
        if (!barnacles[i].equals(line[i + 1])) {
          System.err.println("Error - the two files have different barnacles");
          System.exit(1);
        }
      }
      for (int i = 0; i < anchors.length; i++) {
        if (!anchors[i].equals(line[0])) {
          System.err.println("Error - the two files have different anchors");
          System.exit(1);
        }
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        for (int j = 0; j < barnacles.length; j++) {
          values[1][i * barnacles.length + j] = Double.parseDouble(line[j + 1]);
        }
      }
      reader.close();

      System.out.println(ext.formStr(file1, 30, true) + ext.formStr(file2, 30, true) + " p="
                         + ext.prettyP(new Ttest(values).getPvalue()) + "\t"
                         + ArrayUtils.toStr(Correlation.Pearson(values)));
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + file2 + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + file2 + "\"");
      System.exit(2);
    }

  }

  public static String normalizeDistances(String dir, String distanceFile, double min, double max) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, anchors, barnacles;
    double[] dists;

    if (!new File(dir + ext.rootOf(distanceFile) + "_norm.xln").exists()) {
      try {
        reader = new BufferedReader(new FileReader(dir + distanceFile));
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        anchors = new String[Integer.parseInt(line[0])];
        barnacles = new String[Integer.parseInt(line[1])];
        dists = new double[anchors.length * barnacles.length];
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        for (int i = 0; i < barnacles.length; i++) {
          barnacles[i] = line[i + 1];
        }
        for (int i = 0; i < anchors.length; i++) {
          line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
          anchors[i] = line[0];
          for (int j = 0; j < barnacles.length; j++) {
            dists[i * barnacles.length + j] = Double.parseDouble(line[j + 1]);
          }
        }
        reader.close();

        System.out.print("Transforming data...");
        dists = Transformations.percentileTransform(dists);
        System.out.println("done");

        try {
          writer = Files.openAppropriateWriter(dir + ext.rootOf(distanceFile) + "_norm.xln");
          writer.println(anchors.length + "\t" + barnacles.length);
          writer.println("Anchor\t" + ArrayUtils.toStr(barnacles));
          for (int i = 0; i < anchors.length; i++) {
            writer.print(anchors[i]);
            for (int j = 0; j < barnacles.length; j++) {
              writer.print("\t" + dists[i * barnacles.length + j]);
            }
            writer.println();
          }
          writer.close();
        } catch (Exception e) {
          System.err.println("Error writing normalized distances");
          e.printStackTrace();
        }
      } catch (FileNotFoundException fnfe) {
        fnfe.printStackTrace();
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + dir + distanceFile + "\"");
        ioe.printStackTrace();
        System.exit(2);
      }
    }

    return ext.rootOf(distanceFile) + "_norm.xln";
  }

  public static String matchPairs(String dir, String distanceFile, boolean minMin_not_maxMin,
                                  int iterations) {
    return matchPairs(dir, distanceFile, minMin_not_maxMin, false, iterations);
  }

  /**
   * @param dir directory with file
   * @param distanceFile distance file , computed from {@link MatchSamples#matchMaker}
   * @param minMin_not_maxMin this is a confusing argument
   * @param antiOptimal instead of selecting the optimal match, do the opposite. You most likely
   *          want this to be false
   * @return
   */
  public static String matchPairs(String dir, String distanceFile, boolean minMin_not_maxMin,
                                  boolean antiOptimal, int iterations) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, cases, controls;
    double[][] dists;
    long time;
    double[] mins;
    double[][] finalDists = null;
    int[][] matches = null;
    int iCase, iControl;
    Set<Integer> usedControls = new HashSet<Integer>();
    time = new Date().getTime();
    try {
      reader = new BufferedReader(new FileReader(dir + distanceFile));
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      cases = new String[Integer.parseInt(line[0])];
      controls = new String[Integer.parseInt(line[1])];
      dists = Matrix.doubleMatrix(cases.length, controls.length, -999);
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      for (int i = 0; i < controls.length; i++) {
        controls[i] = line[i + 1];
      }
      for (int i = 0; i < cases.length; i++) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        cases[i] = line[0];
        dists[i] = ArrayUtils.toDoubleArray(ArrayUtils.subArray(line, 1));
      }
      reader.close();
      matches = new int[cases.length][iterations];
      for (int i = 0; i < matches.length; i++) {
        for (int j = 0; j < matches[i].length; j++) {
          matches[i][j] = -1;
        }
      }
      finalDists = new double[cases.length][iterations];
      for (int i = 0; i < finalDists.length; i++) {
        for (int j = 0; j < finalDists[i].length; j++) {
          finalDists[i][j] = -1;
        }
      }
      for (int iter = 0; iter < iterations; iter++) {

        int numDropped = cases.length - controls.length;
        while (findMinTwoD(matches) == -1) {
          // System.out.println(Array.countIf(matches, -1));
          mins = new double[cases.length];
          for (int i = 0; i < cases.length; i++) {
            if (antiOptimal) {
              mins[i] = (matches[i][iter] == -1 ? findMinOrMaxIndex(dists[i], usedControls, false)
                                                : (minMin_not_maxMin ? Double.NEGATIVE_INFINITY
                                                                     : Double.NEGATIVE_INFINITY));
            } else {
              mins[i] = (matches[i][iter] == -1 ? findMinOrMaxIndex(dists[i], usedControls, true)
                                                : (minMin_not_maxMin ? Double.POSITIVE_INFINITY
                                                                     : Double.NEGATIVE_INFINITY));
            }
          }
          if (antiOptimal) {
            iCase = ArrayUtils.maxIndex(mins);
            iControl = ArrayUtils.maxIndex(dists[iCase]);
          } else {
            iCase = minMin_not_maxMin ? ArrayUtils.minIndex(mins) : ArrayUtils.maxIndex(mins);
            if (iCase == -1) break;
            iControl = ArrayUtils.minIndex(dists[iCase]);
            if (iControl == -1) break;
          }
          if ((!minMin_not_maxMin) && (numDropped > 0)) {
            for (int i = 0; i < controls.length; i++) {
              dists[iCase][i] = Double.NEGATIVE_INFINITY;
            }
            System.out.println("Dropped anchor: " + cases[iCase] + ". More cases than controls.");
            numDropped--;
          } else {
            matches[iCase][iter] = iControl;
            finalDists[iCase][iter] = dists[iCase][iControl];
            for (int i = 0; i < cases.length; i++) {
              if (antiOptimal) {
                dists[i][iControl] = Double.NEGATIVE_INFINITY;
              } else {
                dists[i][iControl] = Double.POSITIVE_INFINITY;
              }
            }
          }
          usedControls.add(iControl);
        }
      }
      writer = Files.openAppropriateWriter(dir + ext.rootOf(distanceFile) + "_"
                                           + (minMin_not_maxMin ? "min" : "max") + "Min.xln");
      String header = "Case";
      for (int i = 0; i < iterations; i++) {
        header += "\tControlMatch\tDistance";
      }
      writer.println(header);
      String row;
      for (int i = 0; i < cases.length; i++) {
        row = cases[i];
        for (int j = 0; j < iterations; j++) {
          row += "\t" + controls[matches[i][j]] + "\t" + finalDists[i][j];
        }
        writer.println(row);
      }
      writer.close();
    } catch (FileNotFoundException fnfe) {
      fnfe.printStackTrace();
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + distanceFile + "\"");
      ioe.printStackTrace();
      System.exit(2);
    }
    System.out.println("Created " + ext.rootOf(distanceFile) + "_"
                       + (minMin_not_maxMin ? "min" : "max") + "Min.xln" + " in "
                       + ext.getTimeElapsed(time));

    return ext.rootOf(distanceFile) + "_" + (minMin_not_maxMin ? "min" : "max") + "Min.xln";
  }

  public static double findMinOrMaxIndex(double[] distances, Set<Integer> skipIndices,
                                         boolean findMin) {
    int index = 0;
    double cur;
    if (findMin) {
      double min = Double.POSITIVE_INFINITY;
      for (int i = 0; i < distances.length; i++) {
        if (skipIndices.contains(i)) continue;
        cur = distances[i];
        if (cur < min) {
          min = cur;
          index = i;
        }
      }
      return index;
    } else {
      double max = Double.NEGATIVE_INFINITY;
      for (int i = 0; i < distances.length; i++) {
        if (skipIndices.contains(i)) continue;
        cur = distances[i];
        if (cur > max) {
          max = cur;
          index = i;
        }
      }
      return index;
    }
  }

  public static int findMinTwoD(int[][] in) {
    int min = in[0][0];
    for (int i = 0; i < in.length; i++) {
      for (int j = 0; j < in[i].length; j++) {
        if (in[i][j] < min) {
          min = in[i][j];
        }
      }
    }
    return min;
  }

  public static void evalAgeSex_and_MDS_separately(String dir, String pairings, String refDistances,
                                                   String demofile, String ageHead,
                                                   String genHead) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, anchors, barnacles, refBarns;
    double[] mdsDists, totalDists;
    long time;
    Vector<String> anchs = new Vector<>();
    Vector<String> barns = new Vector<>();
    DoubleVector distV = new DoubleVector();
    int ageIndex, genIndex, numRefAnch, iAnch;
    int[][][] data;
    double[] sumAges, sumGenders;
    String results;
    double mean, stdev;
    int count;

    time = new Date().getTime();

    try {
      reader = new BufferedReader(new FileReader(dir + pairings));
      ext.checkHeader(reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE),
                      new String[] {"Anchor", "BarnaclePair"}, true);
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        anchs.add(line[0]);
        barns.add(line[1]);
        distV.add(Double.parseDouble(line[2]));
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + pairings + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + pairings + "\"");
      System.exit(2);
    }
    anchors = ArrayUtils.toStringArray(anchs);
    barnacles = ArrayUtils.toStringArray(barns);
    totalDists = Doubles.toArray(distV);
    mdsDists = new double[anchors.length];

    try {
      reader = new BufferedReader(new FileReader(dir + refDistances));
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      numRefAnch = Integer.parseInt(line[0]);
      if (anchors.length != Integer.parseInt(line[0])) {
        System.err.println("Warning - number of reference anchors (" + numRefAnch
                           + ") is not the same as the number of anchors (" + anchors.length + ")");
      }
      refBarns = new String[Integer.parseInt(line[1])];
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      for (int i = 0; i < refBarns.length; i++) {
        refBarns[i] = line[i + 1];
      }
      for (int i = 0; i < numRefAnch; i++) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        iAnch = ext.indexOfStr(line[0], anchors);
        mdsDists[iAnch] = Double.parseDouble(line[1 + ext.indexOfStr(barnacles[iAnch], refBarns)]);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + refDistances
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + refDistances + ".xln" + "\"");
      System.exit(2);
    }
    System.out.println("Finished reading in " + refDistances + " in " + ext.getTimeElapsed(time));

    data = new int[anchors.length][2][];
    try {
      reader = new BufferedReader(new FileReader(dir + demofile));
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      ageIndex = ext.indexOfStr(ageHead, line);
      genIndex = ext.indexOfStr(genHead, line);

      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        if (ext.indexOfStr(line[0], anchors) >= 0) {
          data[ext.indexOfStr(line[0], anchors)][0] = new int[] {Integer.parseInt(line[ageIndex]),
                                                                 Integer.parseInt(line[genIndex])};
        } else if (ext.indexOfStr(line[0], barnacles) >= 0) {
          data[ext.indexOfStr(line[0],
                              barnacles)][1] = new int[] {Integer.parseInt(line[ageIndex]),
                                                          Integer.parseInt(line[genIndex])};
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + demofile + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + demofile + "\"");
      System.exit(2);
    }

    sumAges = new double[3];
    sumGenders = new double[3];
    for (int i = 0; i < anchors.length; i++) {
      if (data[i][0] == null) {
        System.err.println("Error - anchor " + anchors[i] + " not found in demofile");
      } else if (data[i][1] == null) {
        System.err.println("Error - barnacle " + barnacles[i] + " not found in demofile");
      }
      sumAges[0] += data[i][0][0];
      sumAges[1] += data[i][1][0];
      sumAges[2] += data[i][0][0] - data[i][1][0];

      sumGenders[0] += data[i][0][1];
      sumGenders[1] += data[i][1][1];
      sumGenders[2] += data[i][0][1] == data[i][1][1] ? 1 : 0;
    }
    for (int i = 0; i < 3; i++) {
      sumAges[i] /= anchors.length;
      sumGenders[i] /= anchors.length;
    }
    System.out.println("Finished evaluating " + pairings + " in " + ext.getTimeElapsed(time));
    System.out.println();

    results = "Cases: " + ext.formPercent(sumGenders[0], 0) + " male, mean AOO="
              + ext.formDeci(sumAges[0], 1, true) + "\n" + "Controls: "
              + ext.formPercent(sumGenders[1], 0) + " male, mean AOO="
              + ext.formDeci(sumAges[1], 1, true) + "\n" + ext.formPercent(sumGenders[2], 0)
              + " gender concordance\n" + "mean age diff b/w case-ctrl: "
              + ext.formDeci(sumAges[2], 1, true) + "\n" + "Mean MDS distance: "
              + ext.formDeci(ArrayUtils.mean(mdsDists), 2) + " (SD: "
              + ext.formDeci(ArrayUtils.stdev(mdsDists), 3, true) + ", range="
              + ext.formDeci(ArrayUtils.min(mdsDists), 2, true) + "-"
              + ext.formDeci(ArrayUtils.max(mdsDists), 2, true) + ")\n";

    mean = ArrayUtils.mean(totalDists);
    stdev = ArrayUtils.stdev(totalDists);
    count = 0;
    sumAges = new double[3];
    sumGenders = new double[3];
    distV = new DoubleVector();
    for (int i = 0; i < anchors.length; i++) {
      if (data[i][0] == null) {
        System.err.println("Error - anchor " + anchors[i] + " not found in demofile");
      } else if (data[i][1] == null) {
        System.err.println("Error - barnacle " + barnacles[i] + " not found in demofile");
      }
      if (totalDists[i] < mean + 3 * stdev) {
        sumAges[0] += data[i][0][0];
        sumAges[1] += data[i][1][0];
        sumAges[2] += data[i][0][0] - data[i][1][0];

        sumGenders[0] += data[i][0][1];
        sumGenders[1] += data[i][1][1];
        sumGenders[2] += data[i][0][1] == data[i][1][1] ? 1 : 0;

        distV.add(mdsDists[i]);

        count++;
      }
    }
    for (int i = 0; i < 3; i++) {
      sumAges[i] /= count;
      sumGenders[i] /= count;
    }
    mdsDists = Doubles.toArray(distV);

    results += "\n\nUsing just the " + count + " of " + anchors.length
               + " pairs that were < 3SD from the mean distance:\n\n";

    results += "Cases: " + ext.formPercent(sumGenders[0], 0) + " male, mean AOO="
               + ext.formDeci(sumAges[0], 1, true) + "\n" + "Controls: "
               + ext.formPercent(sumGenders[1], 0) + " male, mean AOO="
               + ext.formDeci(sumAges[1], 1, true) + "\n" + ext.formPercent(sumGenders[2], 0)
               + " gender concordance\n" + "mean age diff b/w case-ctrl: "
               + ext.formDeci(sumAges[2], 1, true) + "\n" + "Mean MDS distance: "
               + ext.formDeci(ArrayUtils.mean(mdsDists), 2) + " (SD: "
               + ext.formDeci(ArrayUtils.stdev(mdsDists), 3, true) + ", range="
               + ext.formDeci(ArrayUtils.min(mdsDists), 2, true) + "-"
               + ext.formDeci(ArrayUtils.max(mdsDists), 2, true) + ")\n";

    System.out.println(results);

    try {
      writer = Files.openAppropriateWriter(dir + ext.rootOf(pairings) + "_summary1.out");
      writer.println(results);
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + pairings + " summary");
      e.printStackTrace();
    }
  }

  public static void eval(String dir, String pairings, String demofile, String[] factorVars) {
    BufferedReader reader;
    String[] line, anchors, barnacles;
    double[] dists, totalDists, dataList;
    Vector<String> anchs = new Vector<>();
    Vector<String> barns = new Vector<>();
    DoubleVector distV = new DoubleVector();
    double[][][] data;
    double[][] sums;
    String results;
    double mean, stdev;
    int count;
    int[] indices;
    boolean problem;
    boolean[] checkForConcordance;
    String[] factors;

    factors = new String[factorVars.length];
    checkForConcordance = new boolean[factors.length];
    for (int i = 0; i < factors.length; i++) {
      if (factorVars[i].endsWith("=concordance")) {
        factors[i] = factorVars[i].substring(0, factorVars[i].lastIndexOf("="));
        checkForConcordance[i] = true;
      } else {
        factors[i] = factorVars[i];
        checkForConcordance[i] = false;
      }
    }

    try {
      reader = new BufferedReader(new FileReader(dir + pairings));
      ext.checkHeader(reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE),
                      new String[] {"Anchor", "BarnaclePair"}, true);
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        anchs.add(line[0]);
        barns.add(line[1]);
        distV.add(Double.parseDouble(line[2]));
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + pairings + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + pairings + "\"");
      System.exit(2);
    }
    anchors = ArrayUtils.toStringArray(anchs);
    barnacles = ArrayUtils.toStringArray(barns);
    totalDists = Doubles.toArray(distV);

    problem = false;
    data = new double[anchors.length][factors.length][];
    try {
      reader = new BufferedReader(new FileReader(dir + demofile));
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      indices = new int[factors.length];
      for (int i = 0; i < factors.length; i++) {
        indices[i] = ext.indexOfStr(factors[i], line);
        if (indices[i] == -1) {
          System.err.println("Error - could not find factor '" + factors[i]
                             + "' in demographics file '" + demofile + "'");
          problem = true;
        }
      }
      if (problem) {
        System.exit(1);
      }

      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        if (ext.indexOfStr(line[0], anchors) >= 0) {
          dataList = new double[factors.length];
          for (int i = 0; i < indices.length; i++) {
            dataList[i] = Double.parseDouble(line[indices[i]]);
          }
          data[ext.indexOfStr(line[0], anchors)][0] = dataList;
        } else if (ext.indexOfStr(line[0], barnacles) >= 0) {
          dataList = new double[factors.length];
          for (int i = 0; i < indices.length; i++) {
            dataList[i] = Double.parseDouble(line[indices[i]]);
          }
          data[ext.indexOfStr(line[0], barnacles)][1] = dataList;
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + dir + demofile + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + dir + demofile + "\"");
      System.exit(2);
    }

    for (int i = 0; i < anchors.length; i++) {
      if (data[i][0] == null) {
        System.err.println("Error - anchor " + anchors[i] + " not found in demofile");
        problem = true;
      } else if (data[i][1] == null) {
        System.err.println("Error - barnacle " + barnacles[i] + " not found in demofile");
        problem = true;
      }
    }
    if (problem) {
      System.exit(1);
    }

    results = "";

    mean = ArrayUtils.mean(totalDists);
    stdev = ArrayUtils.stdev(totalDists);
    for (int outliers = 0; outliers < 2; outliers++) {
      count = 0;
      distV = new DoubleVector();
      sums = new double[factors.length][3];
      for (int i = 0; i < factors.length; i++) {
        for (int j = 0; j < anchors.length; j++) {
          if (outliers == 0 || totalDists[j] < mean + 2 * stdev) {
            sums[i][0] += data[j][0][i];
            sums[i][1] += data[j][1][i];
            if (checkForConcordance[i]) {
              sums[i][2] += data[j][0][i] == data[j][1][i] ? 1 : 0;
            } else {
              sums[i][2] += data[j][0][i] - data[j][1][i];
            }
            distV.add(totalDists[i]);
            if (i == 0) {
              count++;
            }
          }
        }
      }
      for (int i = 0; i < factors.length; i++) {
        for (int j = 0; j < 3; j++) {
          sums[i][j] /= count;
        }
      }
      dists = Doubles.toArray(distV);

      if (outliers == 0) {
        results += "Using all " + count + " of " + anchors.length + " pairs:\n\n";
      } else {
        results += "\n\nUsing just the " + count + " of " + anchors.length
                   + " pairs that were < 3SD from the mean distance:\n\n";
      }
      results += "\tAnchors\tBarncls\tConcord\tMean Diff\n";
      for (int i = 0; i < factors.length; i++) {
        results += factors[i];
        if (checkForConcordance[i]) {
          results += "\t" + ext.formPercent(sums[i][0], 0) + "\t" + ext.formPercent(sums[i][1], 1)
                     + "\t" + ext.formPercent(sums[i][2], 1) + "\t";
        } else {
          results += "\t" + ext.formDeci(sums[i][0], 1, true) + "\t"
                     + ext.formDeci(sums[i][1], 2, true) + "\t\t"
                     + ext.formDeci(sums[i][2], 2, true);
        }
        results += "\n";
      }
      results += "Mean distance: " + ext.formDeci(ArrayUtils.mean(dists), 2) + " (SD: "
                 + ext.formDeci(ArrayUtils.stdev(dists), 3, true) + ", range="
                 + ext.formDeci(ArrayUtils.min(dists), 2, true) + "-"
                 + ext.formDeci(ArrayUtils.max(dists), 2, true) + ")\n";
    }

    System.out.println(results);
    Files.write(results, dir + ext.rootOf(pairings) + "_summary.out");
  }

  public static void matchFromParameters(String filename, Logger log) {
    List<String> paramV;
    String factorFile, casesFile, controlsFile, demographicsFile, coordsFile;
    String[] line;
    int[] coords;
    String[] factors, demoFactors;
    double[] factorWeights;
    String file, pairs;
    boolean normalize;

    paramV = Files.parseControlFile(filename, "match",
                                    new String[] {"factors.txt normalizeAllFactorsFirst",
                                                  "anchorIDs.txt", "barnacleIDs.txt",
                                                  "demographics.dat Age Sex=concordance",
                                                  "# File with indices for the x-axis and y-axis values to visualize;",
                                                  "factors.txt 1 2",
                                                  "# Factors and weights to use when minimizing distance in C-dimensional space",
                                                  "Age 4", "Sex 1", "PCA1 14", "PCA2 14"},
                                    log);
    if (paramV != null) {
      line = paramV.get(0).split(PSF.Regex.GREEDY_WHITESPACE);
      factorFile = line[0];
      normalize = false;
      for (int i = 1; i < line.length; i++) {
        if (line[i].toLowerCase().startsWith("norm")) {
          normalize = true;
        } else {
          System.err.println("Error - do not know what to do with parameter '" + line[i] + "'");
        }
      }
      casesFile = paramV.get(1);
      controlsFile = paramV.get(2);
      line = paramV.get(3).trim().split(PSF.Regex.GREEDY_WHITESPACE);
      demographicsFile = line[0];
      demoFactors = ArrayUtils.subArray(line, 1);
      line = paramV.get(4).trim().split(PSF.Regex.GREEDY_WHITESPACE);
      coordsFile = line[0];
      coords = new int[] {Integer.parseInt(line[1]), Integer.parseInt(line[2])};
      factors = new String[paramV.size() - 5];
      factorWeights = new double[paramV.size() - 5];
      for (int i = 0; i < factors.length; i++) {
        line = paramV.get(5 + i).trim().split(PSF.Regex.GREEDY_WHITESPACE);
        factors[i] = line[0];
        factorWeights[i] = Double.parseDouble(line[1]);
      }

      // TODO: Update CRF functionality
      log.reportError("CRF Functionality currently disabled for MatchSamples.");
      /*
       * try { file = matchMaker("", casesFile, controlsFile, factorFile, factors, factorWeights,
       * normalize); file = normalizeDistances("", file, 0, 100); pairs = matchPairs("", file,
       * true); eval("", pairs, demographicsFile, demoFactors); new MatchesVisualized("", casesFile,
       * controlsFile, coordsFile, coords, pairs); pairs = matchPairs("", file, false); new
       * MatchesVisualized("", casesFile, controlsFile, coordsFile, coords, pairs); eval("", pairs,
       * demographicsFile, demoFactors); new MatchesVisualized("", casesFile, controlsFile,
       * coordsFile, coords, pairs); ext.waitForResponse(); } catch (Exception e) {
       * log.reportError("Error matching files"); log.reportException(e); }
       */
    }
  }

  public static String buildVisFile(String dir, String fullResultsFile, int currentIteration) {
    String tempVisHelperFile = dir + "vis_helper_" + (currentIteration + 1) + ".temp";
    File f = new File(tempVisHelperFile);
    if (f.exists()) {
      f.delete();
    }
    PrintWriter writer;
    int controlColumn = (currentIteration == 0 ? (currentIteration + 1)
                                               : ((currentIteration * 2) + 1));
    int distanceColumn = controlColumn + 1;
    Vector<String> data = HashVec.loadFileToVec(dir + "/matches/" + fullResultsFile, true, false,
                                                true);
    writer = Files.getAppropriateWriter(tempVisHelperFile);
    String[] line;
    writer.println("Case\tControl\tDistance");
    for (int i = 0; i < data.size(); i++) {
      line = data.elementAt(i).trim().split(PSF.Regex.GREEDY_WHITESPACE);
      writer.println(line[0] + "\t" + line[controlColumn] + "\t" + line[distanceColumn]);
    }
    writer.close();
    return tempVisHelperFile;
  }

  public String run(int iterations, Logger log) {
    String distanceMatrix, pairsFile = null;
    boolean minMin = (this.minOrMax.equalsIgnoreCase("minmin"));

    try {
      String matchesDirectory = this.dir + "matches" + "/";
      Files.ensurePathExists(matchesDirectory);
      distanceMatrix = MatchSamples.matchMaker(matchesDirectory, "/../" + this.samplesFile,
                                               "/../" + this.factorsFile, this.cols, this.loadings,
                                               true);

      String p;
      pairsFile = MatchSamples.matchPairs(matchesDirectory, distanceMatrix, minMin, iterations);

      if (minMin == true) {
        p = pairsFile.split("_minMin")[0] + ".xln";
      } else {
        p = pairsFile.split("_maxMin")[0] + ".xln";
      }
      /*
       * if (!this.skipAge) { MatchSamples.evalAgeSex_and_MDS_separately(matchesDirectory,
       * pairsFile, p, "/../" + factorsFile, "AGE01", "SEX"); MatchSamples.eval(matchesDirectory,
       * pairsFile, "/../" + factorsFile, new String[] {"AGE01", "SEX=concordance"}); }
       */

    } catch (Exception e) {
      e.printStackTrace();
    }
    return pairsFile;

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String d = "./";
    String samplesFile = "samples.txt";
    String factorsFile = "factors.txt";
    String minOrMax = "minmin";
    String[] factorNames = new String[] {"PC1", "PC2"};
    String skipVisualizer = "false";
    String skipAgeSex = "true";
    String resultsToVis = d + "matches/distances_PC1x1,PC2x1_minMin.xln";
    boolean hideExtraControls = false;
    boolean skipVis = false;
    boolean skipAge = true;
    boolean visOnly = false;

    // int[] factorIndices = new int[] {1,2,3,4,5,6,7,8,9,10};
    // int[] factorIndices = new int[] {1,2};

    // String[] factorTargets = new String[] {"C1_norm", "C2_norm",
    // "Age_norm", "AgeAtExam_norm"};
    double[] factorLoadings = new double[] {16, 16, 4, 1};

    String clusterfile = "cluster.genome";
    int iterations = 1;

    String usage = "\n" + "gwas.MatchSamples requires 0-1 arguments\n"
                   + "   (0) directory (i.e. dir=" + d + " (default))\n"
                   + "   (1) sample file with cases (1)  and controls (0) (i.e. samples="
                   + samplesFile + " (default))\n" + "   (3) file with factors (i.e. factors="
                   + factorsFile + " (default))\n"
                   + " (4) column names of factors in clusterfile (i.e. columns="
                   + ArrayUtils.toStr(factorNames, ",") + " (default))\n"
                   + "   (5) clusterfile (i.e. clusterfile=" + clusterfile + " (default))\n"
                   + "   (6) number of control sets to generate (i.e. iterations=1 (default))\n"
                   + "   (7) minMin or maxMin (i.e. minOrMax=minMin (default))\n" + ""
                   + "   (8) skipVisualizer - use this if submitting a non-interactive job (i.e. skipVisualizer=false (default))\n"
                   + "   (9) skipAgeSex - skip matching that uses age and sex (i.e. skipAgeSex=true (default))\n"
                   + "" + "   (10) -visOnly - use this flag to skip matching and only visualize.\n"
                   + ""
                   + "   (11) resultsToVis - used with visOnly flag to point to results file.\n"
                   + ""
                   + "   (12) hideExtraControls - hide unmatched controls in visualizer (i.e. hideExtraControls=false (default)).\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        d = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("samples=")) {
        samplesFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("factors=")) {
        factorsFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("columns=")) {
        factorNames = arg.split("=")[1].split(",");
        numArgs--;
      } else if (arg.startsWith("loadings=")) {
        factorLoadings = ArrayUtils.toDoubleArray(arg.split("=")[1].split(","));
        numArgs--;
      } else if (arg.startsWith("clusterfile=")) {
        clusterfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("iterations=")) {
        iterations = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("minOrMax=")) {
        minOrMax = arg.split("=")[1].toLowerCase();
        numArgs--;
      } else if (arg.startsWith("skipVisualizer=")) {
        skipVisualizer = arg.split("=")[1].toLowerCase();
        numArgs--;
      } else if (arg.startsWith("-visOnly")) {
        visOnly = true;
        numArgs--;
      } else if (arg.startsWith("resultsToVis=")) {
        resultsToVis = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("hideExtraControls=")) {
        hideExtraControls = Boolean.parseBoolean(arg.split("=")[1]);
        numArgs--;
      }
    }

    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    if (visOnly) {
      Logger log = new Logger(d + "vis.log");
      String visFile;
      for (int i = 0; i < iterations; i++) {
        visFile = buildVisFile(d, resultsToVis, i);
        int[] factorIndices = ext.indexFactors(factorNames,
                                               Files.getHeaderOfFile(d + factorsFile,
                                                                     PSF.Regex.GREEDY_WHITESPACE,
                                                                     log),
                                               false, log, true);
        new MatchesVisualized(d, samplesFile, factorsFile,
                              new int[] {factorIndices[0], factorIndices[1]}, visFile,
                              hideExtraControls);

      }
    } else {
      if (factorNames.length < 2) {
        System.err.println("At least two factors are required. Found: "
                           + ArrayUtils.toStr(factorNames, ","));
        System.exit(1);
      }

      if (factorNames.length != factorLoadings.length) {
        System.err.println("Number of loadings and number of factors must match.");
        System.exit(1);
      }

      boolean minMin = true;
      if (minOrMax.equalsIgnoreCase("minmin")) {
        minMin = true;
      } else if (minOrMax.equalsIgnoreCase("maxmin")) {
        minMin = false;
      } else {
        System.err.println(usage);
        System.exit(1);
      }
      if (skipVisualizer.equalsIgnoreCase("true")) {
        skipVis = true;
      }
      if (skipAgeSex.equalsIgnoreCase("false")) {
        skipAge = false;
      }

      MatchSamples match = new MatchSamples(d, samplesFile, factorsFile, factorNames,
                                            factorLoadings, minOrMax, skipVis, skipAge);
      Logger log = new Logger();
      String pairs = match.run(iterations, log);

      if (!skipVis) {
        String visFile;
        for (int i = 0; i < iterations; i++) {
          visFile = buildVisFile(d, pairs, i);
          int[] factorIndices = ext.indexFactors(factorNames,
                                                 Files.getHeaderOfFile(d + factorsFile,
                                                                       PSF.Regex.GREEDY_WHITESPACE,
                                                                       log),
                                                 false, log, true);
          new MatchesVisualized(match.dir, samplesFile, factorsFile,
                                new int[] {factorIndices[0], factorIndices[1]}, visFile,
                                hideExtraControls);

        }
      }
    }

  }
}