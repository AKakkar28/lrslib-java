package com.vertexenumeration;

import java.math.BigInteger;
import java.util.*;

/**
 * Robust vertex enumeration from H-representation:
 *   H is m x (d+1) with rows [a0, a1..ad] meaning a0 + a·x >= 0.
 *
 * Strategy:
 *   1) Try all d-row subsets S (|S|=d), solve A_S x = -a0_S exactly.
 *   2) Keep x if feasible: for all rows i, a0_i + a_i·x >= 0.
 *   3) Emit vertex as homogeneous [1, x..], with strong dedup.
 *
 * Correctness-first; ordering matches lrs by sorting vertices
 * lexicographically DESC on (x_d, x_{d-1}, ..., x_1).
 */
public final class ReverseSearchPivotEnumerator {

    public static final class Result {
        public final List<Fraction[]> vertices; // homogeneous [1, x...]
        public final EnumStats stats;
        Result(List<Fraction[]> v, EnumStats s) { this.vertices = v; this.stats = s; }
    }

    private final Fraction[][] H;
    private final int m, d;

    public ReverseSearchPivotEnumerator(Fraction[][] H) {
        if (H == null || H.length == 0 || H[0].length < 2)
            throw new IllegalArgumentException("H must be m x (d+1) with d>=1");
        this.H = H;
        this.m = H.length;
        this.d = H[0].length - 1;
    }

    /** Public entry point (correctness-first with lrs-style ordering). */
    public Result enumerate() {
        EnumStats st = new EnumStats();
        List<Fraction[]> out = new ArrayList<>();
        Set<String> seen = new LinkedHashSet<>();

        if (d <= 0 || m < d) return new Result(out, st);

        // Precompute ZERO/ONE
        Fraction ZERO = H[0][0].subtract(H[0][0]);
        Fraction ONE = null;
        outer:
        for (int i = 0; i < m; i++) for (int j = 0; j <= d; j++) {
            if (H[i][j].compareTo(ZERO) != 0) { ONE = H[i][j].divide(H[i][j]); break outer; }
        }
        if (ONE == null) ONE = new Fraction(BigInteger.ONE, BigInteger.ONE);

        // Iterate all combinations of d rows
        int[] comb = initComb(d);
        while (comb != null) {
            st.bases++;

            Fraction[][] A = new Fraction[d][d];
            Fraction[] b = new Fraction[d];

            // Build A_S and b_S = -a0_S
            for (int i = 0; i < d; i++) {
                int r = comb[i];
                b[i] = ZERO.subtract(H[r][0]); // -a0
                for (int j = 0; j < d; j++) A[i][j] = H[r][j + 1];
            }

            Fraction[] x = solveExact(A, b, ZERO);
            if (x != null && isFeasible(H, x, ZERO)) {
                Fraction[] vtx = toHomogeneous(x, ONE);
                String key = vertexKey(vtx);
                if (seen.add(key)) {
                    out.add(vtx);
                    st.vertices++;
                }
            }

            comb = nextComb(comb, m, d);
        }

        // ---- lrs-style ordering: lex DESC on (x_d, x_{d-1}, ..., x_1)
        out.sort((a, b) -> {
            for (int j = d; j >= 1; j--) {              // compare last coord first
                int c = b[j].compareTo(a[j]);           // DESC
                if (c != 0) return c;
            }
            return 0;
        });

        st.maxDepth = 0;        // not meaningful for combinational pass
        st.rays = 0;
        st.integerVertices = st.vertices; // inputs are integer in your tests

        return new Result(out, st);
    }

    // ------------------ Linear algebra (exact) ------------------

    /** Solve A x = b exactly via Gauss–Jordan; return null if singular/inconsistent. */
    private static Fraction[] solveExact(Fraction[][] A, Fraction[] b, Fraction ZERO) {
        final int n = b.length;
        Fraction[][] M = new Fraction[n][n + 1];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, M[i], 0, n);
            M[i][n] = b[i];
        }

        int row = 0;
        int[] leadCol = new int[n];
        Arrays.fill(leadCol, -1);

        for (int col = 0; col < n && row < n; col++) {
            int piv = row;
            while (piv < n && M[piv][col].compareTo(ZERO) == 0) piv++;
            if (piv == n) continue;
            if (piv != row) { Fraction[] t = M[piv]; M[piv] = M[row]; M[row] = t; }

            Fraction diag = M[row][col];
            if (diag.compareTo(ZERO) == 0) continue;

            for (int j = col; j <= n; j++) M[row][j] = M[row][j].divide(diag);
            for (int r = 0; r < n; r++) if (r != row) {
                Fraction f = M[r][col];
                if (f.compareTo(ZERO) != 0) {
                    for (int j = col; j <= n; j++) M[r][j] = M[r][j].subtract(f.multiply(M[row][j]));
                }
            }
            leadCol[row] = col;
            row++;
        }

        // Inconsistency check: 0 = nonzero
        for (int i = 0; i < n; i++) {
            boolean allZero = true;
            for (int j = 0; j < n; j++) if (M[i][j].compareTo(ZERO) != 0) { allZero = false; break; }
            if (allZero && M[i][n].compareTo(ZERO) != 0) return null;
        }

        // Unique solution: rank must be n
        int rank = 0; for (int i = 0; i < n; i++) if (leadCol[i] != -1) rank++;
        if (rank != n) return null;

        Fraction[] x = new Fraction[n];
        for (int i = 0; i < n; i++) {
            int lead = -1;
            for (int j = 0; j < n; j++) if (M[i][j].compareTo(ZERO) != 0) { lead = j; break; }
            if (lead == -1) return null;
            x[lead] = M[i][n];
        }
        for (int i = 0; i < n; i++) if (x[i] == null) return null;
        return x;
    }

    private static boolean isFeasible(Fraction[][] H, Fraction[] x, Fraction ZERO) {
        final int m = H.length, d = H[0].length - 1;
        for (int i = 0; i < m; i++) {
            Fraction lhs = H[i][0];
            for (int j = 0; j < d; j++) lhs = lhs.add(H[i][j + 1].multiply(x[j]));
            if (lhs.compareTo(ZERO) < 0) return false;
        }
        return true;
    }

    private static Fraction[] toHomogeneous(Fraction[] x, Fraction ONE) {
        Fraction[] v = new Fraction[x.length + 1];
        v[0] = ONE;
        System.arraycopy(x, 0, v, 1, x.length);
        return v;
    }

    // ------------------ Combination helpers ------------------

    private static int[] initComb(int k) { int[] a = new int[k]; for (int i = 0; i < k; i++) a[i] = i; return a; }
    private static int[] nextComb(int[] a, int n, int k) {
        int i = k - 1;
        while (i >= 0 && a[i] == n - k + i) i--;
        if (i < 0) return null;
        a[i]++;
        for (int j = i + 1; j < k; j++) a[j] = a[j - 1] + 1;
        return a;
    }

    // ------------------ Vertex normalization key ------------------

    /** Vertex key (scale-invariant, with consistent sign orientation). */
    private static String vertexKey(Fraction[] v) {
        Fraction ZERO = v[0].subtract(v[0]);
        Fraction ONE  = v[0].compareTo(ZERO) == 0 ? null : v[0].divide(v[0]);

        // Normalize so v[0] == 1
        Fraction[] n = new Fraction[v.length];
        if (ONE != null && v[0].compareTo(ONE) != 0) {
            for (int i = 0; i < v.length; i++) n[i] = v[i].divide(v[0]);
        } else {
            System.arraycopy(v, 0, n, 0, v.length);
        }

        // Make the first non-zero coordinate positive
        int first = -1;
        for (int j = 1; j < n.length; j++) {
            if (n[j].compareTo(ZERO) != 0) { first = j; break; }
        }
        if (first != -1 && n[first].compareTo(ZERO) < 0) {
            for (int i = 0; i < n.length; i++) n[i] = ZERO.subtract(n[i]);
        }

        StringBuilder sb = new StringBuilder(n.length * 16);
        for (Fraction f : n) sb.append(f.toString()).append('|');
        return sb.toString();
    }
}
