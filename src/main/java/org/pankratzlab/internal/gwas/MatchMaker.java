package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.internal.gwas.FactorLoadings;
import org.pankratzlab.kdmatch.KDMatch;
import org.pankratzlab.kdmatch.KDTree;
import org.pankratzlab.kdmatch.Match;
import org.pankratzlab.kdmatch.Sample;
import org.pankratzlab.kdmatch.SelectOptimizedNeighbors;

import com.google.common.collect.Sets;

public class MatchMaker {

	private static List<Match> kdMatchMaker(Path baseDir, List<Sample> caseList, List<Sample> controlList,
			int initialNumSelect, int finalNumSelect, FactorLoadings factorLoadings, int threads, Logger log) {

		KDTree<Sample> kdTree = new KDTree<>(factorLoadings.getDoubleLoadings().size());
		log.info("Assuming 1 ID column and " + (factorLoadings.getDoubleLoadings().size()) + " data columns");
		log.info("Building tree for: " + caseList.get(1).getGroup());

		KDTree.addSamplesToTree(kdTree, controlList.stream());

		log.info("selecting initial " + initialNumSelect + " nearest neighbors for " + caseList.get(1).getGroup());

		List<Match> naiveMatches = KDTree.getNearestNeighborsForSamples(kdTree, caseList.stream(), initialNumSelect)
				.collect(Collectors.toList());

		String outputBase = baseDir + File.separator + caseList.get(1).getGroup() + "match.naive.txt.gz";
		log.info("reporting full baseline selection of " + initialNumSelect + " nearest neighbors to " + outputBase);
		try {
			KDMatch.writeToFile(naiveMatches.stream(), outputBase,
					factorLoadings.getDoubleFactorNames().toArray(new String[1]),
					factorLoadings.getDoubleFactorNames().toArray(new String[1]), initialNumSelect);
		} catch (IOException e) {
			e.printStackTrace();
		}

		String outputOpt = baseDir + File.separator + caseList.get(1).getGroup() + "match.optimized.txt.gz";

		log.info("selecting optimized nearest neighbors");

		List<Match> optimizedMatches = null;
		try {
			optimizedMatches = SelectOptimizedNeighbors.optimizeDuplicates(naiveMatches, finalNumSelect,
					threads, log).collect(Collectors.toList());
			log.info("reporting optimized selection of " + finalNumSelect + " nearest neighbors to " + outputOpt);
			KDMatch.writeToFile(optimizedMatches.stream(), outputOpt, factorLoadings.getDoubleFactorNames().toArray(new String[1]),
					factorLoadings.getDoubleFactorNames().toArray(new String[1]), finalNumSelect);
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		} catch (ExecutionException e2) {
			e2.printStackTrace();
		} catch (IOException e3) {
			e3.printStackTrace();
		}
		
		return optimizedMatches;

	}

	public static void runMatching(Path dir, Path inputSamples, FactorLoadings factorLoadings, int initialNumSelect,
			int finalNumSelect, int threads, Logger log) throws IOException {
		int[] numericColumnsToUseForClustering = getLoadingIndices(inputSamples, factorLoadings.getFactors(), false,
				log);
		int[] factorColumnsToAssignGroup = getLoadingIndices(inputSamples, factorLoadings.getFactors(), true, log);
		int idColumn = 0;
		Map<String, List<Sample>> casesGroupedByStringFactor = new HashMap<String, List<Sample>>();
		Map<String, List<Sample>> controlsGroupedByStringFactor = new HashMap<String, List<Sample>>();
		Files.lines(inputSamples).map(l -> l.split("\t")).skip(1).map(k -> Sample.parseSample(k, idColumn,
				numericColumnsToUseForClustering, factorColumnsToAssignGroup, factorLoadings))
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

		List<Match> results = casesGroupedByStringFactor.entrySet().stream()
				.map(e -> kdMatchMaker(dir, e.getValue(), controlsGroupedByStringFactor.get(e.getKey()),
						initialNumSelect, finalNumSelect, factorLoadings, threads, log))
				.flatMap(List::stream).collect(Collectors.toList());
	}

	private static int[] getLoadingIndices(Path sampleFile, Map<String, String> factorLoadings, boolean force,
			Logger log) {
		Set<String> forcedCols = new LinkedHashSet<String>();
		Set<String> doubleCols = new LinkedHashSet<String>();
		double temp;
		for (Entry<String, String> s : factorLoadings.entrySet()) {
			if (FactorLoadings.ACCEPTED_LOADINGS.contains(s.getValue())) {
				forcedCols.add(s.getKey());
			} else {
				try {
					temp = Double.parseDouble(s.getValue());
					doubleCols.add(s.getKey());
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
		try (BufferedReader r = org.pankratzlab.common.Files
				.getAppropriateReader(sampleFile.toAbsolutePath().toString())) {
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

	public static void main(String[] args) {

		Path d = Paths.get("./");
		FactorLoadings factorLoadings = null;
		Path samples = Paths.get(d + "samples.txt");
		boolean normalize = true;
		int numArgs = args.length;
		int finalNumSelect = 4;
		int threads = Runtime.getRuntime().availableProcessors();
		Logger log;

		String usage = "\n" + "gwas.MatchMaker requires at least 1 argument\n"
				+ "(1) Working directory (e.g. dir=./ (default))\n"
				+ "(2) Samples filename with ID column, case/control column (1/0), and factor columns (e.g. samples=samples.txt (default))\n"
				+ "(3) Factors columns to use and their loadings; *keep in file order (e.g. factors=PC1:4,PC2:4,sex:force,phenograph:force)\n"
				+ "(4) Normalize inputs (e.g. normalize=true (default))\n"
				+ "(5) Iterations - number of controls to match to each case (e.g. iterations=4 (default))\n";

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
				factorLoadings = new FactorLoadings(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("normalize=")) {
				normalize = Boolean.parseBoolean(arg.split("=")[1]);
			} else if (arg.startsWith("iterations=")) {
				finalNumSelect = Integer.parseInt(arg.split("=")[1]);
			}
		}
		int initialNumSelect = finalNumSelect * 5;
		samples = Paths.get(d.toString() + "\\" + samples.toString());
		log = Logger.getAnonymousLogger();
		log.info("Starting sample match using k-d tree nearest neighbors.");

		try {
			runMatching(d, samples, factorLoadings, initialNumSelect, finalNumSelect, threads, log);
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
