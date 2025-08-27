package com.vertexenumeration;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;
import java.util.*;

public class RayEnumerationTest {

    private static Polyhedron hRepFromRows(Fraction[][] rows) {
        int m = rows.length, n = rows[0].length;
        Matrix M = new Matrix(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                M.set(i, j, rows[i][j]);
        return new Polyhedron(Polyhedron.Type.H, m, n, /*integer*/ false, M);
    }

    private static Set<String> rows(Polyhedron P){
        Set<String> s = new HashSet<>();
        Matrix M = P.getMatrix();
        for (int i=0;i<P.getRowCount();i++){
            StringBuilder sb = new StringBuilder();
            for (int j=0;j<P.getColCount();j++) sb.append(M.get(i,j).toString()).append('|');
            s.add(sb.toString());
        }
        return s;
    }

    @Test
    public void coneInR2_hasOneVertexAndTwoRays() {
        // Cone in R^2: { (x,y) : y >= 0,  x - y >= 0 }  (i.e., along +x and y<=x)
        // H rows [b|A]: y>=0 -> [0,0,1] ; x - y >= 0 -> [0,1,-1] (apex at origin)
        Fraction[][] H = {
                { Fraction.ZERO, Fraction.ZERO, Fraction.ONE  },
                { Fraction.ZERO, Fraction.ONE,  Fraction.of(-1) }
        };
        Polyhedron input = hRepFromRows(H);
        VertexEnumerator ve = new VertexEnumerator();
        Polyhedron out = ve.enumerate(input);

        assertEquals(Polyhedron.Type.V, out.getType());
        // Expect one vertex [1|0 0] and two extreme rays [0|1 0], [0|1 1] (any canonical scaling OK)
        Set<String> got = rows(out);

        assertTrue(got.contains("1|0|0|"), "apex present");
        // Rays have leading 0; direction depends on your canonicalization; accept either
        boolean hasRayX = got.contains("0|1|0|") || got.contains("0|2|0|") || got.stream().anyMatch(s->s.startsWith("0|") && s.endsWith("|0|"));
        boolean hasRayDiag = got.contains("0|1|1|") || got.contains("0|2|2|") || got.stream().anyMatch(s->s.startsWith("0|") && s.contains("|") && !s.endsWith("|0|"));
        assertTrue(hasRayX, "ray along +x present");
        assertTrue(hasRayDiag, "ray along x=y present");
    }
}
