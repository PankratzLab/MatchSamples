package org.pankratzlab.internal.gwas;

import org.pankratzlab.common.ext;

public class MatchingVariable {
  final public String prettyName;
  final public String headerName;
  final public boolean isBinary;

  private int headerIndex = -1;

  private ContinuousDataBox continuousDataBox;
  private BinaryDataBox binaryDataBox;

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

  public void setContinuousDataBox(ContinuousDataBox box) {
    cantBeBinary();
    this.continuousDataBox = box;
  }

  public void setBinaryDataBox(BinaryDataBox box) {
    mustBeBinary();
    this.binaryDataBox = box;
  }

  public double getCaseAvg() {
    cantBeBinary();
    return continuousDataBox.getCaseAvg(this);
  }

  public double getControlAvg() {
    cantBeBinary();
    return continuousDataBox.getControlAvg(this);
  }

  public double getConcordance() {
    mustBeBinary();
    return binaryDataBox.getConcordance(this);
  }

  private void cantBeBinary() {
    if (isBinary) {
      throw new UnsupportedOperationException("This operation is not supported for binary matching variables");
    }
  }

  private void mustBeBinary() {
    if (!isBinary) {
      throw new UnsupportedOperationException("This operation is not supported for non-binary matching variables");
    }
  }
}
