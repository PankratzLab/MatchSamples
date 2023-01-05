package org.pankratzlab.internal.gwas;


import org.junit.Test;
import org.junit.jupiter.api.function.Executable;

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

  @Test
  public void testMethodRestrictionBasedOnBinary() {
    MatchingVariable binary = new MatchingVariable("foo", true);
    Executable[] nonBinaryOnlyMethods = new Executable[] {() -> binary.addToControlSum(0.2),
                                                          () -> binary.addToCaseSum(0.3),
                                                          binary::getControlAvg,
                                                          binary::getCaseAvg};
    for (Executable method : nonBinaryOnlyMethods) {
      assertThrows(UnsupportedOperationException.class, method);
    }

    MatchingVariable continuous = new MatchingVariable("bar", false);
    Executable[] binaryOnlyMethods = new Executable[] {() -> continuous.recordControlValue("case",
                                                                                           1),
                                                       () -> continuous.setCaseValue("case", 0),
                                                       continuous::getConcordance};
    for (Executable method : binaryOnlyMethods) {
      assertThrows(UnsupportedOperationException.class, method);
    }
  }

}
