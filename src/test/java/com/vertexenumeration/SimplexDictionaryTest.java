package com.vertexenumeration;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;
import java.util.*;

public class SimplexDictionaryTest {

    /** Unit square in R^2:
     *  H rows are b + A x >= 0 for:
     *   x >= 0, y >= 0, 1 - x >= 0, 1 - y >= 0
     *   → rows: [0,1,0], [0,0,1], [1,-1,0], [1,0,-1]
     */
    private static Fraction[][] unitSquareH() {
        Fraction Z = Fraction.ZERO, O = Fraction.ONE;
        return new Fraction[][]{
                { Z,  O,  Z },
                { Z,  Z,  O },
                { O,  Fraction.of(-1),  Z },
                { O,  Z,  Fraction.of(-1) }
        };
    }

    @Test
    public void testVertexAndSlacks() {
        Fraction[][] H = unitSquareH();
        // Basis rows for vertex (0,0): inequalities 0:x>=0 and 1:y>=0 tight
        int[] basis = {0,1};
        SimplexDictionary D = new SimplexDictionary(H, basis);
        Fraction[] v = D.vertex();
        assertEquals(Fraction.ZERO, v[0]);
        assertEquals(Fraction.ZERO, v[1]);
        // Slacks: for 1-x>=0 and 1-y>=0 should be 1 at (0,0)
        assertEquals(Fraction.ONE, D.slack(2));
        assertEquals(Fraction.ONE, D.slack(3));
    }

    @Test
    public void testChildrenBasesAndParent() {
        Fraction[][] H = unitSquareH();
        // Start at (0,0) as above
        SimplexDictionary root = new SimplexDictionary(H, new int[]{0,1});
        List<int[]> kids = root.childrenBases();
        // Neighbors should be bases for (1,0) and (0,1):
        //   replace row0 (x>=0) with row2 (1-x>=0)  → (1,0)
        //   replace row1 (y>=0) with row3 (1-y>=0)  → (0,1)
        Set<String> got = new HashSet<>();
        for (int[] b : kids) got.add(Arrays.toString(b));
        assertTrue(got.contains(Arrays.toString(new int[]{1,2}))    // {1,2} sorted
                || got.contains(Arrays.toString(new int[]{2,1})));
        assertTrue(got.contains(Arrays.toString(new int[]{0,3}))
                || got.contains(Arrays.toString(new int[]{3,0})));

        // Check parent mapping is consistent (root is parent of its neighbors)
        for (int[] childB : kids) {
            SimplexDictionary cd = new SimplexDictionary(H, childB);
            int[] par = cd.parentBasis();
            assertNotNull(par);
            assertArrayEquals(new int[]{0,1}, par);
        }
    }

    @Test
    public void testLexOrderOfChildrenIsDeterministic() {
        Fraction[][] H = unitSquareH();
        SimplexDictionary root = new SimplexDictionary(H, new int[]{0,1});
        List<int[]> kids = root.childrenBases();
        // Should be lex-sorted ascending
        int[] prev = null;
        for (int[] b : kids) {
            if (prev != null) {
                for (int i = 0; i < b.length; i++) {
                    int c = Integer.compare(prev[i], b[i]);
                    if (c < 0) break;
                    if (c > 0) fail("childrenBases not in lex order");
                }
            }
            prev = b;
        }
    }
}
