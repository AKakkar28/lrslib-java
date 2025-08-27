package com.vertexenumeration;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;
import java.util.*;
import java.io.*;
import java.nio.charset.StandardCharsets;

public class VertexEnumerationEndToEndTest {

    private static Polyhedron hRepFromRows(Fraction[][] rows) {
        int m = rows.length, n = rows[0].length;
        Matrix M = new Matrix(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                M.set(i, j, rows[i][j]);
        return new Polyhedron(Polyhedron.Type.H, m, n, /*integer*/ false, M);
    }

    private static Set<String> asRowSet(Polyhedron vRep) {
        Set<String> s = new HashSet<>();
        Matrix M = vRep.getMatrix();
        for (int i = 0; i < vRep.getRowCount(); i++) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < vRep.getColCount(); j++) {
                sb.append(M.get(i, j).toString()).append('|');
            }
            s.add(sb.toString());
        }
        return s;
    }

    @Test
    public void square_H_to_V() {
        // Unit square in R^2: x≥0, y≥0, 1−x≥0, 1−y≥0
        Fraction Z = Fraction.ZERO, O = Fraction.ONE;
        Fraction[][] H = {
                { Z,  O,  Z },             // x >= 0
                { Z,  Z,  O },             // y >= 0
                { O,  Fraction.of(-1), Z}, // 1 - x >= 0
                { O,  Z, Fraction.of(-1)}  // 1 - y >= 0
        };

        Polyhedron input = hRepFromRows(H);
        VertexEnumerator ve = new VertexEnumerator();
        Polyhedron out = ve.enumerate(input);
        assertEquals(Polyhedron.Type.V, out.getType());
        assertEquals(3, out.getColCount()); // [lead | x y]
        // Expected 4 vertices, no rays
        assertTrue(out.getRowCount() >= 4);
        Set<String> got = asRowSet(out);

        Set<String> want = new HashSet<>(Arrays.asList(
                "1|0|0|", "1|1|0|", "1|0|1|", "1|1|1|"
        ));
        assertTrue(got.containsAll(want), "All square vertices present");
    }

    @Test
    public void cube_H_to_V() {
        // Unit cube in R^3: x,y,z in [0,1]
        Fraction Z = Fraction.ZERO, O = Fraction.ONE;
        Fraction[][] H = {
                { Z,  O,  Z,  Z },               // x >= 0
                { Z,  Z,  O,  Z },               // y >= 0
                { Z,  Z,  Z,  O },               // z >= 0
                { O,  Fraction.of(-1), Z,  Z },  // 1 - x >= 0
                { O,  Z,  Fraction.of(-1), Z },  // 1 - y >= 0
                { O,  Z,  Z,  Fraction.of(-1) }  // 1 - z >= 0
        };

        Polyhedron input = hRepFromRows(H);
        VertexEnumerator ve = new VertexEnumerator();
        Polyhedron out = ve.enumerate(input);
        assertEquals(Polyhedron.Type.V, out.getType());
        assertEquals(4, out.getColCount()); // [lead | x y z]
        assertTrue(out.getRowCount() >= 8);

        // all 8 {0,1}^3 vertices
        Set<String> want = new HashSet<>();
        for (int x=0;x<=1;x++)
            for (int y=0;y<=1;y++)
                for (int z=0;z<=1;z++)
                    want.add("1|" + x + "|" + y + "|" + z + "|");

        Set<String> got = asRowSet(out);
        assertTrue(got.containsAll(want), "All cube vertices present");
    }
}
