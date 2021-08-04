package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.StringJoiner;

import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.PSF;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

public class Eval {
  public LinkedHashMap<String, String> evalArgs;
  public ArrayList<String> quantFactorNames;
  public Path results;
  public int iterations;

  public Eval(LinkedHashMap<String, String> evalArgs, ArrayList<String> quantFactorNames,
              Path results, int iterations) {
    this.evalArgs = evalArgs;
    this.quantFactorNames = quantFactorNames;
    this.results = results;
    this.iterations = iterations;
  }

  public void evaluateFinalResults() {
    try (BufferedReader resultsReader = Files.getAppropriateReader(results.toString())) {
      PrintWriter writer = Files.getAppropriateWriter(results.getRoot().toString()
                                                      + "evalResults.txt");
      Multimap<Integer, Integer> relatedQuantIndices = HashMultimap.create();
      Multimap<Integer, Integer> relatedEvalIndices = HashMultimap.create();
      String[] header = resultsReader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      for (int i = 0; i < quantFactorNames.size(); i++) {
        for (int j = 0; j < header.length; j++) {
          if (header[j].contains(quantFactorNames.get(i))) {
            relatedQuantIndices.put(i, j);
          }
        }
      }
      ArrayList<String> evalNames = new ArrayList<String>();
      for (String s : evalArgs.keySet()) {
        evalNames.add(s);
      }
      for (int i = 0; i < evalNames.size(); i++) {
        for (int j = 0; j < header.length; j++) {
          if (header[j].contains(evalNames.get(i))) {
            relatedEvalIndices.put(i, j);
          }
        }
      }

      String line = resultsReader.readLine();
      while (line != null) {
        String[] lineArr = line.trim().split(PSF.Regex.GREEDY_WHITESPACE);
        for (Integer x : relatedQuantIndices.keySet()) {
          for (Integer y : relatedQuantIndices.get(x)) {

          }
        }
      }

      StringJoiner joiner = new StringJoiner("\t");
      for (String s : quantFactorNames) {
        for (int i = 1; i <= iterations; i++) {
          joiner.add(s + "_round_" + i);
        }
      }
      for (String s : evalNames) {
        for (int i = 1; i <= iterations; i++) {
          joiner.add(s + "_round_" + i);
        }
      }

    } catch (IOException ioe) {
      System.out.println("Results file " + results.toString() + " not found.");
    }

  }

}
