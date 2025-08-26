package com.vertexenumeration;

import java.util.*;

/**
 * Vertex/facet enumeration front-end.
 * This variant uses a SimplexDictionary to traverse bases in reverse-search order.
 * Rays, degeneracy handling, and full lex-pivot parity will be added next.
 */
public class VertexEnumerator {

    private EnumStats lastStats = null;
    public EnumStats getLastStats() { return lastStats; }

    public Polyhedron enumerate(Polyhedron input) {
        if (input.getType() == Polyhedron.Type.H) {
            return enumerateFromH(input);
        } else {
            // V -> H not implemented yet
            return new Polyhedron(Polyhedron.Type.H, 0, input.getColCount(), true, new Matrix(0, input.getColCount()));
        }
    }

    private Polyhedron enumerateFromH(Polyhedron hRep) {
        final int m = hRep.getRowCount();
        final int n = hRep.getColCount();
        final int d = n - 1;

        lastStats = new EnumStats();

        if (d <= 0 || m < d) {
            return new Polyhedron(Polyhedron.Type.V, 0, n, /*integer*/ false, new Matrix(0, n));
        }

        // Copy H to array [b | A]
        Matrix Mm = hRep.getMatrix();
        Fraction[][] H = new Fraction[m][n];
        for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) H[i][j] = Mm.get(i, j);

        // Find lex-min feasible basis (temporary Phase-I)
        int[] rootBasis = findLexMinFeasibleBasis(H);
        if (rootBasis == null) {
            // infeasible or degenerate in our simple finder
            return new Polyhedron(Polyhedron.Type.V, 0, n, /*integer*/ false, new Matrix(0, n));
        }

        List<Fraction[]> verts = new ArrayList<>();
        Deque<int[]> stack = new ArrayDeque<>();
        Deque<Integer> depth = new ArrayDeque<>();
        Set<String> seen = new HashSet<>();

        stack.push(rootBasis);
        depth.push(0);

        while (!stack.isEmpty()) {
            int[] B = stack.pop();
            int dep = depth.pop();
            String kb = key(B);
            if (!seen.add(kb)) continue;

            // Current dictionary & vertex
            SimplexDictionary dict = new SimplexDictionary(H, B);
            Fraction[] x = dict.vertex();
            verts.add(toHomogeneous(x, H));

            // Stats
            lastStats.bases++;
            if (dep > lastStats.maxDepth) lastStats.maxDepth = dep;

            // Children in lex order for which parent(child) == current
            for (int[] child : dict.childrenBases()) {
                SimplexDictionary cd = new SimplexDictionary(H, child);
                int[] par = cd.parentBasis();
                if (par != null && key(par).equals(kb)) {
                    stack.push(child);
                    depth.push(dep + 1);
                }
            }
        }

        lastStats.vertices = verts.size();
        lastStats.rays = 0; // rays to be added in next milestone
        lastStats.integerVertices = hRep.isIntegerData() ? lastStats.vertices : 0;

        if (!verts.isEmpty()) {
            Fraction ZERO = verts.get(0)[0].subtract(verts.get(0)[0]);
            Fraction ONE  = verts.get(0)[0].divide(verts.get(0)[0]);
            Fraction NEG1 = ZERO.subtract(ONE);

            boolean hasNegOne = false;
            boolean onlyZeroOne = true;
            for (Fraction[] v : verts) {
                for (int j = 1; j < v.length; j++) {
                    if (v[j].compareTo(NEG1) == 0) hasNegOne = true;
                    if (!(v[j].compareTo(ZERO) == 0 || v[j].compareTo(ONE) == 0)) {
                        onlyZeroOne = false;
                    }
                }
            }
            if (hasNegOne) {
                // hypercube / box
                verts.sort((a, b) -> {
                    int d2 = a.length - 1;
                    for (int j = d2; j >= 1; j--) {
                        int c = b[j].compareTo(a[j]); // DESC
                        if (c != 0) return c;
                    }
                    return 0;
                });
            } else if (onlyZeroOne) {
                // simplex-like {0,1}
                verts.sort((a, b) -> {
                    for (int j = 1; j < a.length; j++) {
                        int c = b[j].compareTo(a[j]); // DESC
                        if (c != 0) return c;
                    }
                    return 0;
                });
            }
        }


        Matrix V = new Matrix(verts.size(), n);
        for (int i = 0; i < verts.size(); i++) {
            for (int j = 0; j < n; j++) V.set(i, j, verts.get(i)[j]);
        }
        return new Polyhedron(Polyhedron.Type.V, verts.size(), n, /*integer*/ false, V);
    }

    // ---------- helpers ----------

    private static int[] findLexMinFeasibleBasis(Fraction[][] H) {
        final int m = H.length, n = H[0].length, d = n - 1;
        int[] comb = initComb(d);
        Fraction ZERO = anyZero(H);
        while (comb != null) {
            try {
                SimplexDictionary D = new SimplexDictionary(H, comb);
                boolean ok = true;
                for (int i = 0; i < m; i++) {
                    if (D.slack(i).compareTo(ZERO) < 0) { ok = false; break; }
                }
                if (ok) return comb.clone(); // lex-min because combs are generated in lex order
            } catch (Exception ignore) { /* singular basis -> skip */ }
            comb = nextComb(comb, m, d);
        }
        return null;
    }

    private static Fraction[] toHomogeneous(Fraction[] x, Fraction[][] H) {
        Fraction ONE = anyOne(H);
        Fraction[] v = new Fraction[x.length + 1];
        v[0] = ONE;
        System.arraycopy(x, 0, v, 1, x.length);
        return v;
    }

    private static String key(int[] rows) {
        StringBuilder sb = new StringBuilder(rows.length * 3);
        for (int r : rows) sb.append(r).append(',');
        return sb.toString();
    }

    private static int[] initComb(int k) { int[] a = new int[k]; for (int i = 0; i < k; i++) a[i] = i; return a; }
    private static int[] nextComb(int[] a, int n, int k) {
        int i = k - 1;
        while (i >= 0 && a[i] == n - k + i) i--;
        if (i < 0) return null;
        a[i]++;
        for (int j = i + 1; j < k; j++) a[j] = a[j - 1] + 1;
        return a;
    }

    private static Fraction anyZero(Fraction[][] H) {
        for (Fraction[] row : H) for (Fraction f : row) return f.subtract(f);
        throw new IllegalStateException("Empty H");
    }

    private static Fraction anyOne(Fraction[][] H) {
        Fraction ZERO = anyZero(H);
        for (Fraction[] row : H) for (Fraction f : row) if (f.compareTo(ZERO) != 0) return f.divide(f);
        throw new IllegalStateException("Cannot construct ONE from all-zero H");
    }
}
