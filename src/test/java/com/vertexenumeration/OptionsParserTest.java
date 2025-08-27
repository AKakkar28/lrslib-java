// src/test/java/com/vertexenumeration/OptionsParserTest.java
package com.vertexenumeration;
import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

public class OptionsParserTest {
    @Test
    public void parsesBasicFlags() {
        String[] args = {
                "-v", "-redund", "-minrep", "-printcobasis", "-seed", "42",
                "-maxdepth","100", "-threads","4", "-eliminate","2,3", "-project","1,2",
                "-linearity","0,5", "square.ine"
        };
        OptionsParser.Parsed p = OptionsParser.parse(args);
        LrsDat d = p.dat;
        assertEquals("square.ine", p.inputPath);
        assertEquals(LrsDat.Mode.VE, d.mode);
        assertTrue(d.redund); assertTrue(d.minrep); assertTrue(d.printCobasis);
        assertEquals(42L, d.seed); assertEquals(100, d.maxDepth); assertEquals(4, d.threads);
        assertArrayEquals(new int[]{1,2}, d.project);  // parser sorts? (ours preserves order—adjust if needed)
        assertArrayEquals(new int[]{2,3}, d.eliminate);
        assertArrayEquals(new int[]{0,5}, d.linearities);
    }
}
