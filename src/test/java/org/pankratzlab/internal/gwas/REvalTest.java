package org.pankratzlab.internal.gwas;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class REvalTest {

  @Test
  public void basicFunctionTest() throws IOException {
    MatchingVariable foo = new MatchingVariable("foo", false);
    MatchingVariable bar = new MatchingVariable("bar", true);
    File status = new File("src/test/resources/status.tsv");
    File phenotype = new File("src/test/resources/phenotype.tsv");

    REval rEval = new REval(new MatchingVariable[] {foo, bar}, status, phenotype);
    rEval.readPhenotypeFile();

    assertEquals(-0.5, foo.getCaseAvg());
    assertEquals(0.5, foo.getControlAvg());
    assertEquals(0.5, bar.getConcordance());
  }
}
