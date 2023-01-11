package org.pankratzlab.internal.gwas;

import org.pankratzlab.common.ext;

public class MatchingVariable {
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
}
