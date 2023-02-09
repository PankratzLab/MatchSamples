package org.pankratzlab.internal.gwas;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.pankratzlab.common.ext;

public class MatchingVariable {
  public static final String NA = "NA";
  public final String prettyName;
  public final String headerName;

  private boolean isBinary = true;

  private int headerIndex = -1;

  private DataBox dataBox;

  private Set<Double> uniqueValues = new HashSet<>();
  private double maxValue = Double.NEGATIVE_INFINITY;
  private double minValue = Double.POSITIVE_INFINITY;

  public MatchingVariable(String headerName) {
    this.headerName = headerName;
    this.prettyName = headerName;
  }

  public int findIndexInHeader(String[] header) {
    headerIndex = ext.indexOfStr(headerName, header);
    return headerIndex;
  }

  public void checkBinary(Double value) {
    if (isBinary) {
      uniqueValues.add(value);
      if (uniqueValues.size() > 2) {
        isBinary = false;
        return;
      }
      maxValue = Math.max(maxValue, value);
      minValue = Math.min(minValue, value);
      if (maxValue - minValue > 1) {
        isBinary = false;
      }
    }
  }

  public boolean isBinary() {
    return isBinary;
  }

  public boolean isContinuous() {
    return !isBinary;
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

  public static MatchingVariable[] fromCommaSeparatedString(String commaSep) {
    String[] names = commaSep.split(",");
    return Arrays.stream(names).map(String::strip).map(MatchingVariable::new)
                 .toArray(MatchingVariable[]::new);
  }

  public static MatchingVariable[] fromNames(Collection<String> names) {
    return names.stream().map(MatchingVariable::new).toArray(MatchingVariable[]::new);
  }
}
