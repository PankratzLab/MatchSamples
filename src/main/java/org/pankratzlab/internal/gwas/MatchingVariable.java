package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;

public class MatchingVariable {
  public static final String NA = "NA";
  final public String prettyName;
  final public String headerName;
  final public boolean isBinary;

  private int headerIndex = -1;

  private DataBox dataBox;

  public MatchingVariable(String headerName, boolean isBinary) {
    this.headerName = headerName;
    this.prettyName = headerName;
    this.isBinary = isBinary;
  }

  public int findIndexInHeader(String[] header) {
    headerIndex = ext.indexOfStr(headerName, header);
    return headerIndex;
  }

  public int getHeaderIndex() {
    return headerIndex;
  }

  public void setDataBox(DataBox dataBox) {
    this.dataBox = dataBox;
  }

  public double getCaseAvg() {
    cantBeBinary();
    return dataBox.getCaseAvg(this);
  }

  public double getControlAvg() {
    cantBeBinary();
    return dataBox.getControlAvg(this);
  }

  public double getConcordance() {
    mustBeBinary();
    return dataBox.getConcordance(this);
  }

  public double getUnivariateP() {
    return dataBox.getUnivariateP(this);
  }

  public double getMultivariateP() {
    return dataBox.getMultivariateP(this);
  }

  public String getTableLine() {
    String caseAvg = this.isBinary ? NA : prettyDecimal(this.getCaseAvg());
    String controlAvg = this.isBinary ? NA : prettyDecimal(this.getControlAvg());
    String concordance = this.isBinary ? prettyDecimal(this.getConcordance()) : NA;
    String univariateP = prettyDecimal(this.getUnivariateP());
    String multivariateP = prettyDecimal(this.getMultivariateP());
    return String.join("\t", this.headerName, caseAvg, controlAvg, concordance, univariateP,
                       multivariateP);
  }

  private static String prettyDecimal(double d) {
    return ext.formDeci(d, 5);
  }

  public void cantBeBinary() {
    if (isBinary) {
      throw new UnsupportedOperationException("This operation is not supported for binary matching variables");
    }
  }

  public void mustBeBinary() {
    if (!isBinary) {
      throw new UnsupportedOperationException("This operation is not supported for non-binary matching variables");
    }
  }

  private static MatchingVariable fromString(String s) {
    String[] values = s.strip().split("\t");
    if (values.length != 2) {
      throw new IllegalArgumentException("Matching variable data appears to have the wrong number of values.");
    }
    String name = values[0].strip();
    boolean binary = Boolean.parseBoolean(values[1].trim());
    return new MatchingVariable(name, binary);
  }

  public static MatchingVariable[] fromFile(File file) throws IOException {
    if (!file.isFile()) {
      throw new IllegalArgumentException("Given matching variable file does not exist or is not a normal file.");
    }
    BufferedReader reader = Files.getAppropriateReader(file.toString());
    String expectedHeader = "header_name\tis_binary";
    String actualHeader = reader.readLine().strip();
    if (!actualHeader.equals(expectedHeader)) {
      throw new IllegalArgumentException("The provided matching variable file has an incorrect header.");
    }
    List<MatchingVariable> matchingVariables = new ArrayList<>();
    String inputLine = reader.readLine();
    while (inputLine != null) {
      matchingVariables.add(fromString(inputLine));
      inputLine = reader.readLine();
    }
    return matchingVariables.toArray(MatchingVariable[]::new);
  }
}
