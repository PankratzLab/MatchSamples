package org.pankratzlab.internal.gwas;

import java.util.Map;

import org.junit.Test;
import org.junit.jupiter.api.function.Executable;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

public class MatchingVariableAndDataBoxTest {

  @Test
  public void testRecordControlValues() {
    MatchingVariable mv = new MatchingVariable("foo", true);
    mv.findIndexInHeader(new String[] {"id", "foo"});
    Map<String, String> pairings = Map.of("cont1", "case1", "cont2", "case1");
    DataBox dataBox = new DataBox(new MatchingVariable[] {}, new MatchingVariable[] {mv}, pairings);

    dataBox.recordData(new String[] {"cont1", "1"});
    dataBox.recordData(new String[] {"cont2", "0"});
    dataBox.recordData(new String[] {"case1", "0"});

    assertEquals(0.5, mv.getConcordance());
  }

  @Test
  public void testCaseValueNotSet() {
    MatchingVariable mv = new MatchingVariable("foo", true);
    mv.findIndexInHeader(new String[] {"id", "foo"});
    Map<String, String> pairings = Map.of("cont1", "case1", "cont2", "case1");

    DataBox dataBox = new DataBox(new MatchingVariable[] {}, new MatchingVariable[] {mv}, pairings);

    dataBox.recordData(new String[] {"cont1", "1"});
    dataBox.recordData(new String[] {"cont2", "0"});

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
