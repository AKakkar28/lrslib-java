package com.vertexenumeration;

import java.util.*;

/**
 * Reverse search enumerator (Java translation of lrslib/reverse.c).
 * Enumerates all vertices of an H-polyhedron via dictionary pivoting.
 */
final class ReverseSearchEnumerator {

    static final class Result {
        final List<Fraction[]> vertices;
        final EnumStats stats;
        Result(List<Fraction[]> v, EnumStats s) { vertices = v; stats = s; }
    }

    private final Fraction[][] H; // inequalities: rows [b | A]
    private final int m, n, d;

    ReverseSearchEnumerator(Fraction[][] H) {
        this.H = H;
        this.m = H.length;
        this.n = (m == 0 ? 0 : H[0].length);
        this.d = Math.max(0, n - 1);
        if (m == 0 || n == 0 || d == 0)
            throw new IllegalArgumentException("Empty/degenerate matrix");
    }

    /** Run reverse search enumeration. */
    Result run() {
        EnumStats stats = new EnumStats();
        List<Fraction[]> verts = new ArrayList<>();

        // Step 1: find lexicographic minimum feasible basis
        int[] rootBasis = findLexMinFeasibleBasis(H);
        if (rootBasis == null)
            return new Result(Collections.emptyList(), stats);

        // Step 2: DFS traversal
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

            SimplexDictionary dict = new SimplexDictionary(H, B);
            Fraction[] x = dict.vertex();
            Fraction[] homog = toHomogeneous(x);
            verts.add(homog);

            stats.bases++;
            stats.vertices++;
            if (isIntegerVertex(homog)) stats.integerVertices++;
            if (dep > stats.maxDepth) stats.maxDepth = dep;

            for (int[] child : dict.childrenBases()) {
                SimplexDictionary cd = new SimplexDictionary(H, child);
                int[] par = cd.parentBasis();
                if (par != null && key(par).equals(kb)) {
                    stack.push(child);
                    depth.push(dep + 1);
                }
            }
        }

        verts.sort((a, b) -> {
            // lrs prints vertices in descending order by coordinates,
            // starting from the last coordinate toward the first
            for (int j = a.length - 1; j >= 1; j--) {
                int cmp = b[j].compareTo(a[j]); // notice: reversed
                if (cmp != 0) return cmp;
            }
            return 0;
        });

        return new Result(verts, stats);
    }

    // ---------- Utilities ----------

    /** Homogenize x → [1, x...] */
    private static Fraction[] toHomogeneous(Fraction[] x) {
        Fraction[] v = new Fraction[x.length + 1];
        v[0] = Fraction.ONE;
        System.arraycopy(x, 0, v, 1, x.length);
        return v;
    }

    /** Check if homogeneous vertex is integer. */
    private static boolean isIntegerVertex(Fraction[] v) {
        for (int j = 1; j < v.length; j++) {
            if (!v[j].denominator().equals(java.math.BigInteger.ONE))
                return false;
        }
        return true;
    }

    /** Lexicographically minimal feasible basis search (brute force). */
    private static int[] findLexMinFeasibleBasis(Fraction[][] H) {
        final int m = H.length, n = H[0].length, d = n - 1;
        int[] comb = initComb(d);
        while (comb != null) {
            try {
                SimplexDictionary D = new SimplexDictionary(H, comb);
                boolean ok = true;
                for (int i = 0; i < m; i++) {
                    if (D.slack(i).compareTo(Fraction.ZERO) < 0) {
                        ok = false;
                        break;
                    }
                }
                if (ok) return comb.clone();
            } catch (Exception ignore) {}
            comb = nextComb(comb, m, d);
        }
        return null;
    }

    private static String key(int[] rows) {
        StringBuilder sb = new StringBuilder(rows.length * 3);
        for (int r : rows) sb.append(r).append(',');
        return sb.toString();
    }

    private static int[] initComb(int k) {
        int[] a = new int[k];
        for (int i = 0; i < k; i++) a[i] = i;
        return a;
    }

    private static int[] nextComb(int[] a, int n, int k) {
        int i = k - 1;
        while (i >= 0 && a[i] == n - k + i) i--;
        if (i < 0) return null;
        a[i]++;
        for (int j = i + 1; j < k; j++) a[j] = a[j - 1] + 1;
        return a;
    }
}
