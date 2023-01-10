package org.pankratzlab.internal.gwas;

import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.junit.jupiter.api.function.Executable;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

public class MatchingVariableAndBinaryDataBoxTest {

  @Test
  public void testRecordControlValues() {
    MatchingVariable mv = new MatchingVariable("foo", true);
    Map<String, String> pairings = Map.of("cont1", "case1", "cont2", "case1");
    BinaryDataBox binaryDataBox = new BinaryDataBox(new MatchingVariable[] {mv}, 3, pairings,
                                                    Set.of("case1"));

    binaryDataBox.recordValue("cont1", mv, 1);
    binaryDataBox.recordValue("cont2", mv, 0);
    binaryDataBox.recordValue("case1", mv, 0);

    assertEquals(0.5, mv.getConcordance());
  }

  @Test
  public void testExceptionIfCaseValueNotSet() {
    MatchingVariable mv = new MatchingVariable("foo", true);
    Map<String, String> pairings = Map.of("cont1", "case1", "cont2", "case1");
    BinaryDataBox binaryDataBox = new BinaryDataBox(new MatchingVariable[] {mv}, 3, pairings,
                                                    Set.of("case1"));

    binaryDataBox.recordValue("cont1", mv, 1);
    binaryDataBox.recordValue("cont2", mv, 0);

    assertThrows(IllegalStateException.class, mv::getConcordance);
  }

  @Test
  public void testMethodRestrictionBasedOnBinary() {
    MatchingVariable binary = new MatchingVariable("foo", true);
    Executable[] nonBinaryOnlyMethods = new Executable[] {binary::getControlAvg,
                                                          binary::getCaseAvg};
    for (Executable method : nonBinaryOnlyMethods) {
      assertThrows(UnsupportedOperationException.class, method);
    }

    MatchingVariable continuous = new MatchingVariable("bar", false);
    Executable[] binaryOnlyMethods = new Executable[] {continuous::getConcordance};
    for (Executable method : binaryOnlyMethods) {
      assertThrows(UnsupportedOperationException.class, method);
    }
  }

}
