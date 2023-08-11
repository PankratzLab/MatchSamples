package org.pankratzlab.internal.gwas;

import java.util.Map;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.function.Executable;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class MatchingVariableAndDataBoxTest {

  @Test
  public void testRecordControlValues() {
    MatchingVariable mv = new MatchingVariable("foo");
    mv.findIndexInHeader(new String[] {"id", "foo"});
    Map<String, String> pairings = Map.of("cont1", "case1", "cont2", "case1");
    DataBox dataBox = new DataBox(new MatchingVariable[] {mv}, pairings);

    dataBox.recordData(new String[] {"cont1", "1"});
    dataBox.recordData(new String[] {"cont2", "0"});
    dataBox.recordData(new String[] {"case1", "0"});

    assertEquals(0.5, mv.getConcordance());
  }

  @Test
  public void testCaseValueNotSet() {
    MatchingVariable mv = new MatchingVariable("foo");
    mv.findIndexInHeader(new String[] {"id", "foo"});
    Map<String, String> pairings = Map.of("cont1", "case1", "cont2", "case1");

    DataBox dataBox = new DataBox(new MatchingVariable[] {mv}, pairings);

    dataBox.recordData(new String[] {"cont1", "1"});
    dataBox.recordData(new String[] {"cont2", "0"});

    assertThrows(IllegalStateException.class, mv::getConcordance);
  }

  @Test
  public void testMethodRestrictionBasedOnBinary() {
    MatchingVariable binary = new MatchingVariable("foo");
    MatchingVariable continuous = new MatchingVariable("bar");

    String[] header = new String[] {"id", "foo", "bar"};
    binary.findIndexInHeader(header);
    continuous.findIndexInHeader(header);

    Map<String, String> pairings = Map.of("cont1", "case1", "cont2", "case1");

    DataBox dataBox = new DataBox(new MatchingVariable[] {binary, continuous}, pairings);

    // MVs are binary by default so this should work with no data recorded
    Executable[] nonBinaryOnlyMethods = new Executable[] {binary::getControlAvg,
                                                          binary::getCaseAvg};
    for (Executable method : nonBinaryOnlyMethods) {
      assertThrows(UnsupportedOperationException.class, method);
    }

    // MVs are binary by default so we have to record a bit of data to make one continuous
    dataBox.recordData(new String[] {"cont1", "0", "2.5"});
    dataBox.recordData(new String[] {"cont2", "1", "1.3"});
    dataBox.recordData(new String[] {"case1", "1", "-3.4"});

    Executable[] binaryOnlyMethods = new Executable[] {continuous::getConcordance};
    for (Executable method : binaryOnlyMethods) {
      assertThrows(UnsupportedOperationException.class, method);
    }
  }

  @Test
  public void testBinaryDetectionBasedOnUniqueValueCount() {
    MatchingVariable mv = new MatchingVariable("foo");
    mv.findIndexInHeader(new String[] {"id", "foo"});

    Map<String, String> pairings = Map.of("cont1", "case1", "cont2", "case1", "cont3", "case1");
    DataBox dataBox = new DataBox(new MatchingVariable[] {mv}, pairings);

    dataBox.recordData(new String[] {"cont1", "1"});
    dataBox.recordData(new String[] {"cont2", "0"});
    dataBox.recordData(new String[] {"case1", "0"});

    assertTrue(mv.isBinary());

    dataBox.recordData(new String[] {"cont3", "1.1"});

    assertFalse(mv.isBinary());
  }

  @Test
  public void testBinaryDetectionBasedOnValueRange() {
    MatchingVariable mv = new MatchingVariable("foo");
    mv.findIndexInHeader(new String[] {"id", "foo"});

    Map<String, String> pairings = Map.of("cont1", "case1", "cont2", "case1", "cont3", "case1");
    DataBox dataBox = new DataBox(new MatchingVariable[] {mv}, pairings);

    dataBox.recordData(new String[] {"cont1", "0"});
    dataBox.recordData(new String[] {"cont2", "0"});
    dataBox.recordData(new String[] {"case1", "0"});

    assertTrue(mv.isBinary());

    dataBox.recordData(new String[] {"cont3", "1.1"});

    assertFalse(mv.isBinary());
  }
}
