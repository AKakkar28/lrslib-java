package com.vertexenumeration;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;
import java.util.*;

public class FacetEnumerationEndToEndTest {

    private static Polyhedron vRepFromRows(Fraction[][] rows) {
        int m = rows.length, n = rows[0].length;
        Matrix M = new Matrix(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                M.set(i, j, rows[i][j]);
        return new Polyhedron(Polyhedron.Type.V, m, n, /*integer*/ false, M);
    }

    private static Set<String> asRowSet(Polyhedron hRep) {
        Set<String> s = new HashSet<>();
        Matrix M = hRep.getMatrix();
        for (int i = 0; i < hRep.getRowCount(); i++) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < hRep.getColCount(); j++) {
                sb.append(M.get(i, j).toString()).append('|');
            }
            s.add(sb.toString());
        }
        return s;
    }

    @Test
    public void square_V_to_H_facets() {
        // vertices of unit square in V-form: [1|x y]
        Fraction[][] V = {
                { Fraction.ONE, Fraction.ZERO, Fraction.ZERO },
                { Fraction.ONE, Fraction.ONE,  Fraction.ZERO },
                { Fraction.ONE, Fraction.ZERO, Fraction.ONE  },
                { Fraction.ONE, Fraction.ONE,  Fraction.ONE  }
        };
        Polyhedron input = vRepFromRows(V);

        Polyhedron out = FacetEnumerator.fromV(input);
        assertEquals(Polyhedron.Type.H, out.getType());
        assertEquals(3, out.getColCount()); // [a0 a1 a2]

        // Expected supporting halfspaces (up to positive scaling):
        //  x >= 0      -> [0, 1, 0]
        //  y >= 0      -> [0, 0, 1]
        //  1 - x >= 0  -> [1, -1, 0]
        //  1 - y >= 0  -> [1, 0, -1]
        Set<String> want = new HashSet<>(Arrays.asList(
                "0|1|0|", "0|0|1|", "1|-1|0|", "1|0|-1|"
        ));

        Set<String> got = asRowSet(out);
        // Accept canonical positive scalings only (FacetEnumerator normalizes first nonzero to +1)
        assertTrue(got.containsAll(want), "All square facets present");
        assertEquals(4, out.getRowCount());
    }
}
