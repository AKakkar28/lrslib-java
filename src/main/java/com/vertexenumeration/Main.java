package com.vertexenumeration;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

public class Main {
    public static void main(String[] args) {
        if (args.length != 1) {
            System.err.println("Usage: java -cp target/<jar or classes> com.vertexenumeration.Main <input-file>");
            System.exit(1);
        }
        String filename = args[0];
        long t0 = System.nanoTime();

        try {
            Polyhedron input = Polyhedron.readFromFile(filename);

            if (input.getType() == Polyhedron.Type.H) {
                ReverseSearchPivotEnumerator algo = new ReverseSearchPivotEnumerator(extractH(input));
                ReverseSearchPivotEnumerator.Result res = algo.enumerate();

                Polyhedron vrep = toV(res.vertices, input.getColCount());
                try (PrintWriter out = new PrintWriter(System.out)) { vrep.write(out); }
                System.out.printf("*Totals: vertices=%d rays=%d bases=%d%n", res.stats.vertices, 0, res.stats.bases);

            } else {
                Polyhedron hrep = FacetEnumerator.fromV(input);
                try (PrintWriter out = new PrintWriter(System.out)) { hrep.write(out); }
            }

            double secs = (System.nanoTime() - t0) / 1_000_000_000.0;
            System.out.printf("*Time=%.3fs%n", secs);

        } catch (FileNotFoundException e) {
            System.err.println("File not found: " + filename);
            System.exit(1);
        } catch (IOException e) {
            System.err.println("I/O error: " + e.getMessage());
            System.exit(1);
        }
    }

    private static Fraction[][] extractH(Polyhedron h) {
        int m = h.getRowCount(), n = h.getColCount();
        Fraction[][] H = new Fraction[m][n];
        Matrix M = h.getMatrix();
        for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) H[i][j] = M.get(i, j);
        return H;
    }

    private static Polyhedron toV(java.util.List<Fraction[]> verts, int ncols) {
        Matrix out = new Matrix(verts.size(), ncols);
        for (int i = 0; i < verts.size(); i++) {
            Fraction[] row = verts.get(i);
            for (int j = 0; j < ncols; j++) out.set(i, j, row[j]);
        }
        return new Polyhedron(Polyhedron.Type.V, verts.size(), ncols, false, out);
    }
}
