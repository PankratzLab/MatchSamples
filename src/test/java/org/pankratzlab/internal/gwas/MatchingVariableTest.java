package org.pankratzlab.internal.gwas;

import org.junit.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

public class MatchingVariableTest {

  @Test
  public void testRecordControlValues() {
    MatchingVariable mv = new MatchingVariable("foo", true);
    String caseId = "case1234";
    mv.recordControlValue(caseId, 1);
    mv.recordControlValue(caseId, 0);

    mv.setCaseValue(caseId, 0);

    assertEquals(0.5, mv.getConcordance());
  }

  @Test
  public void testExceptionIfCaseValueNotSet() {
    MatchingVariable mv = new MatchingVariable("foo", true);
    String caseId = "case1234";
    mv.recordControlValue(caseId, 1);
    mv.recordControlValue(caseId, 0);

    assertThrows(IllegalStateException.class, mv::getConcordance);
  }

  @Test
  public void testControlsRecordedAfterCaseSet() {
    MatchingVariable mv = new MatchingVariable("foo", true);
    String caseId = "case1234";
    mv.recordControlValue(caseId, 1);
    mv.setCaseValue(caseId, 0);
    mv.recordControlValue(caseId, 0);

    assertEquals(0.5, mv.getConcordance());
  }
}
