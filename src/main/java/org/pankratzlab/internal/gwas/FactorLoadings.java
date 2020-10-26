package org.pankratzlab.internal.gwas;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

public class FactorLoadings {
	public LinkedHashMap<String, String> factorLoadings;
	
	public static final Set<String> ACCEPTED_LOADINGS = new HashSet<>(Arrays.asList("force", "binary", "nominal"));
	
	public FactorLoadings(String factorsArgString) {
		factorLoadings = readArgString(factorsArgString);
	}
	
	public LinkedHashMap<String, String> getFactors(){
		return this.factorLoadings;
	}
	
	public ArrayList<Double> getDoubleLoadings(){
		ArrayList<Double> doubleLoadings = new ArrayList<Double>();
		for (String s : factorLoadings.values()) {
			if (ACCEPTED_LOADINGS.contains(s)) {
				continue;
			} else {
				try {
					doubleLoadings.add(Double.parseDouble(s));
				} catch (NumberFormatException nfe) {
					nfe.printStackTrace();
					System.exit(1);
				}
			}
		}
		return doubleLoadings;
	}
	
	public ArrayList<String> getDoubleFactorNames(){
		ArrayList<String> doubleFactorNames = new ArrayList<String>();
		for (String s : factorLoadings.keySet()) {
			if (ACCEPTED_LOADINGS.contains(s)) {
				continue;
			} else {
				doubleFactorNames.add(s);
			}
		}
		return doubleFactorNames;
	}
	public static LinkedHashMap<String, String> readArgString(String factorsArgString) {
		LinkedHashMap<String, String> fl = new LinkedHashMap<String, String>();
		String[] commaSplit = factorsArgString.split(",");
		for (String factorLoading : commaSplit) {
			if (factorLoading.contains(":")) {
				String[] colonSplit = factorLoading.split(":");
				if (!ACCEPTED_LOADINGS.contains(colonSplit[1].toLowerCase())) {
					try {
						fl.put(colonSplit[0], Double.toString(Double.parseDouble(colonSplit[1])));
					} catch(NumberFormatException nfe) {
						System.err.println("Loading " + colonSplit[1] + " not recognized.");
					}
				} else {
					fl.put(colonSplit[0], colonSplit[1]);
				}
			}
		}
		return fl;
	}
}