package org.pankratzlab.internal.gwas;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

public class FactorLoadings {
  public LinkedHashMap<String, String> factorLoadings;

  public static final Set<String> ACCEPTED_LOADINGS = new HashSet<>(Arrays.asList("force", "binary",
                                                                                  "nominal"));

  public FactorLoadings(String factorsArgString) {
    factorLoadings = readArgString(factorsArgString);
  }

  public LinkedHashMap<String, String> getFactors() {
    return this.factorLoadings;
  }

  public ArrayList<Double> getDoubleLoadings() {
    ArrayList<Double> doubleLoadings = new ArrayList<Double>();
    for (String s : factorLoadings.values()) {
      if (ACCEPTED_LOADINGS.contains(s)) {
        continue;
      } else {
        try {
          if (s.contains("nom_")) {
            doubleLoadings.add(Double.parseDouble(s.trim().split("_")[1]));
          } else {
            doubleLoadings.add(Double.parseDouble(s));
          }
        } catch (NumberFormatException nfe) {
          nfe.printStackTrace();
          System.exit(1);
        }
      }
    }
    return doubleLoadings;
  }

  public ArrayList<String> getQuantFactorNames() {
    ArrayList<String> quantFactorNames = new ArrayList<String>();
    for (Entry<String, String> s : factorLoadings.entrySet()) {
      if (ACCEPTED_LOADINGS.contains(s.getValue())) {
        continue;
      } else {
        if (!s.getValue().contains("nom_")) {
          quantFactorNames.add(s.getKey());
        }
      }
    }
    return quantFactorNames;
  }

  public ArrayList<String> getNumericFactorNames() {
    ArrayList<String> numericFactorNames = new ArrayList<String>();
    for (Entry<String, String> s : factorLoadings.entrySet()) {
      if (ACCEPTED_LOADINGS.contains(s.getValue())) {
        continue;
      } else {
        numericFactorNames.add(s.getKey());
      }
    }
    return numericFactorNames;
  }

  public ArrayList<String> getNumericFactorNames(boolean includeNominal) {
    ArrayList<String> numericFactorNames = new ArrayList<String>();
    for (Entry<String, String> s : factorLoadings.entrySet()) {
      if (ACCEPTED_LOADINGS.contains(s.getValue())) {
        continue;
      } else {
        if (s.getValue().contains("nom_")) {
          if (includeNominal) {
            numericFactorNames.add(s.getKey());
          }
        } else {
          numericFactorNames.add(s.getKey());
        }
      }
    }
    return numericFactorNames;
  }

  public ArrayList<String> getNominalFactorNames() {
    ArrayList<String> nominalFactorNames = new ArrayList<String>();
    for (Entry<String, String> s : factorLoadings.entrySet()) {
      if (s.getValue().contains("nom_")) {
        nominalFactorNames.add(s.getKey());
      }
    }
    return nominalFactorNames;
  }

  public static LinkedHashMap<String, String> readArgString(String factorsArgString) {
    LinkedHashMap<String, String> fl = new LinkedHashMap<String, String>();
    String[] commaSplit = factorsArgString.split(",");
    double defaultWeight = 1;
    int index = 0;
    for (String factorLoading : commaSplit) {
      if (factorLoading.contains(":")) {
        String[] colonSplit = factorLoading.split(":");
        if (colonSplit.length == 3 || colonSplit[1].equalsIgnoreCase("nominal")) {
          if (colonSplit[1].equalsIgnoreCase("nominal") && !(colonSplit.length == 3)) {
            fl.put(colonSplit[0], "nom_" + Double.toString(defaultWeight));
          } else if (colonSplit.length == 3 && colonSplit[1].equalsIgnoreCase("nominal")) {
            try {
              fl.put(colonSplit[0], "nom_" + Double.toString(Double.parseDouble(colonSplit[2])));
            } catch (NumberFormatException nfe) {
              System.err.println("Nominal loading " + colonSplit[2] + " not recognized.");
              System.exit(1);
            }
          } else {
            try {
              fl.put(colonSplit[0], "nom_" + Double.toString(Double.parseDouble(colonSplit[1])));
            } catch (NumberFormatException nfe) {
              System.err.println("Nominal loading " + colonSplit[1] + " not recognized.");
              System.exit(1);
            }
          }
        } else {
          if (!ACCEPTED_LOADINGS.contains(colonSplit[1].toLowerCase())) {
            try {
              fl.put(colonSplit[0], Double.toString(Double.parseDouble(colonSplit[1])));
            } catch (NumberFormatException nfe) {
              System.err.println("Loading " + colonSplit[1] + " not recognized.");
              System.exit(1);
            }
          } else {
            fl.put(colonSplit[0], colonSplit[1]);
          }
        }
      }
      index++;
    }
    return fl;
  }
}