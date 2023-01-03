package org.pankratzlab.internal.gwas;

import org.pankratzlab.common.ext;

public class MatchingVariable {
  final String prettyName;
  final String headerName;
  final boolean isBinary;

  private int headerIndex = -1;

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
}
