package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.HashMap;

import org.pankratzlab.common.Files;

import com.google.common.collect.HashBasedTable;

public class MatchSamplesPhenoFiles {

  public static void createFiles(String dir, String sampleFile, String criteriaFile) {

    BufferedReader sampleReader;
    BufferedReader criteriaReader;
    PrintWriter writer;
    HashBasedTable<Integer, String, String> sampleInfo = HashBasedTable.create();
    HashMap<Integer, String> indexMap = new HashMap<Integer, String>();

    try {
      sampleReader = new BufferedReader(new FileReader(dir + sampleFile));

      String[] header = sampleReader.readLine().trim().split("\t");
      for (int i = 0; i < header.length; i++) {
        indexMap.put(i, header[i]);
      }

      int row = 0;
      String line = sampleReader.readLine();
      while (line != null) {
        String[] entries = line.trim().split("\t");
        for (int i = 0; i < entries.length; i++) {
          sampleInfo.put(row, indexMap.get(i), entries[i]);
        }
        row++;
        line = sampleReader.readLine();
      }

      criteriaReader = new BufferedReader(new FileReader(dir + criteriaFile));
      String critLine = criteriaReader.readLine();
      int numPCs = 4;
      while (critLine != null) {
        HashMap<String, String> criteriaMap = new HashMap<String, String>();
        String[] critLineArray = critLine.trim().split("\t");
        String filename = "";
        boolean sexCol = true;
        boolean ageCol = false;
        for (int i = 0; i < critLineArray.length; i++) {
          String column = critLineArray[i].split(":")[0];
          String filters = critLineArray[i].split(":")[1];
          if (column.equals("ImputedRace")) {
            if (filters.equals("Hispanic") || filters.equals("Asian")) {
              numPCs = 3;
            } else if (filters.equals("White")) {
              numPCs = 2;
            } else {
              numPCs = 4;
            }
          }
          if (column.equals("location")) {
            if (filters.contains("ovary") || filters.contains("testis")) {
              sexCol = false;
            }
          }
          if (column.equals("matched_sex")) {
            sexCol = false;
          }
          if (column.equals("age_cat")) {
            ageCol = true;
          }
          criteriaMap.put(column, filters);
          for (String filter : criteriaMap.get(column).split(",")) {
            if (filename.length() != 0)
              filename += ("_" + filter);
            else
              filename += filter;
          }
        }
        if (!sexCol) {
          filename += "_sex.pheno";
        } else {
          if (ageCol) {
            filename += "_ageGrp";
          }
          filename += ".pheno";
        }
        writer = Files.openAppropriateWriter(dir + filename);
        boolean hit = true;
        int numHits = 0;
        boolean multipleStudies = false;
        for (Integer currRow : sampleInfo.rowKeySet()) {
          for (String filterColumn : criteriaMap.keySet()) {
            String filterList = criteriaMap.get(filterColumn).trim();
            if (filterColumn.equals("SoftMatch") && filterList.split(",").length > 1) {
              multipleStudies = true;
            }
            if (!filterList.contains(sampleInfo.get(currRow, filterColumn))) {
              hit = false;
            }
          }
          if (hit) {
            if (numHits == 0) {
              writer.write("fid" + "\t" + "iid" + "\t" + "case_control");
              if (sexCol) {
                writer.write("\t" + "sex");
              }
              for (int i = 1; i <= numPCs; i++) {
                writer.write("\t" + "pc" + i);
              }
              if (multipleStudies) {
                writer.write("\t" + "PG2_bin_soft" + "\t" + "Michigan_bin_soft" + "\t"
                             + "California_bin_soft" + "\n");
              } else {
                writer.write("\n");
              }
            }
            writer.write(sampleInfo.get(currRow, "fid") + "\t" + sampleInfo.get(currRow, "iid")
                         + "\t" + sampleInfo.get(currRow, "case_control"));
            if (sexCol) {
              writer.write("\t" + sampleInfo.get(currRow, "sex"));
            }
            for (int i = 1; i <= numPCs; i++) {
              writer.write("\t" + sampleInfo.get(currRow, "pc" + i));
            }
            if (multipleStudies) {
              writer.write("\t" + sampleInfo.get(currRow, "PG2_bin_soft") + "\t"
                           + sampleInfo.get(currRow, "Michigan_bin_soft") + "\t"
                           + sampleInfo.get(currRow, "California_bin_soft") + "\n");
            } else {
              writer.write("\n");
            }
            numHits++;
          }
          hit = true;
        }
        writer.close();
        critLine = criteriaReader.readLine();
      }

    } catch (Exception e) {
      e.printStackTrace();
      System.exit(1);
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;

    String usage = "\\n" + "gwas.MatchSamplesPhenoFiles requires 3 arguments\n"
                   + "(1) directory where files are located (dir= )\n"
                   + "(2) file with samples to include relative to dir (samples= )\n"
                   + "(3) file with tab-delimited criteria for multiple pheno files (criteria= )\n"
                   + "";

    String dir = "";
    String sampleFile = "";
    String criteriaFile = "";
    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("samples=")) {
        sampleFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("criteria=")) {
        criteriaFile = arg.split("=")[1];
        numArgs--;
      }
    }

    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    createFiles(dir, sampleFile, criteriaFile);

  }
}