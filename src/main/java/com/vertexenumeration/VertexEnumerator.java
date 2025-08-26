package com.vertexenumeration;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Minimal exact vertex enumerator for H-representations.
 * For dimension d = n-1, try all m choose d subsets of tight constraints,
 * solve A x = -a0, and keep feasible solutions (a0 + a·x >= 0 for all rows).
 * Produces a V-representation (vertices only; rays not handled yet).
 *
 * This is a correctness baseline so cube.ine prints 8 vertices.
 * Later, replace with reverse-search (Avis–Fukuda) for scale and parity.
 */
public class VertexEnumerator {

    public Polyhedron enumerate(Polyhedron input) {
        if (input.getType() == Polyhedron.Type.H) {
            return enumerateFromH(input);
        } else {
            // V -> H not implemented yet
            Matrix empty = new Matrix(0, input.getColCount());
            return new Polyhedron(Polyhedron.Type.H, 0, input.getColCount(), true, empty);
        }
    }

    private Polyhedron enumerateFromH(Polyhedron hRep) {
        final int m = hRep.getRowCount();   // number of inequalities
        final int n = hRep.getColCount();   // 1 + dimension
        final int d = n - 1;

        if (d <= 0 || m < d) {
            Matrix empty = new Matrix(0, n);
            return new Polyhedron(Polyhedron.Type.V, 0, n, true, empty);
        }

        Matrix M = hRep.getMatrix();
        Fraction[][] H = toArray(M, m, n);

        // Build reusable ZERO/ONE without relying on Fraction.zero()/one()
        final Fraction ZERO = anyZero(H);
        final Fraction ONE  = anyOne(H);

        List<Fraction[]> verts = new ArrayList<>();
        Set<String> seen = new HashSet<>();

        // Iterate all d-combinations of row indices
        int[] comb = initComb(d);
        while (comb != null) {
            // A x = b, where for each chosen row r: a0 + a·x = 0  =>  a·x = -a0
            Fraction[][] A = new Fraction[d][d];
            Fraction[] b = new Fraction[d];

            for (int i = 0; i < d; i++) {
                int r = comb[i];
                // -a0  (we don't assume negate(); compute 0 - a0)
                b[i] = ZERO.subtract(H[r][0]);
                for (int j = 0; j < d; j++) {
                    A[i][j] = H[r][j + 1];
                }
            }

            Fraction[] x = solve(A, b, ZERO);
            if (x != null && feasible(H, x, ZERO)) {
                // homogeneous vertex: [1, x1, ..., xd]
                Fraction[] v = new Fraction[n];
                v[0] = ONE; // 1
                for (int j = 0; j < d; j++) v[j + 1] = x[j];

                String key = canonical(v);
                if (seen.add(key)) verts.add(v);
            }

            comb = nextComb(comb, m, d);
        }

        Matrix V = new Matrix(verts.size(), n);
        for (int i = 0; i < verts.size(); i++) {
            for (int j = 0; j < n; j++) {
                V.set(i, j, verts.get(i)[j]);
            }
        }
        return new Polyhedron(Polyhedron.Type.V, verts.size(), n, /* integerData */ false, V);

        // If you track integer flag on P: use hRep.isInteger() instead of 'true'.
    }

    // ---------- helpers ----------

    private static Fraction[][] toArray(Matrix M, int rows, int cols) {
        Fraction[][] out = new Fraction[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                out[i][j] = M.get(i, j);
            }
        }
        return out;
    }

    /** Return a ZERO of the same Fraction type using x - x. */
    private static Fraction anyZero(Fraction[][] H) {
        for (Fraction[] row : H) {
            for (Fraction f : row) {
                // f - f = 0
                return f.subtract(f);
            }
        }
        throw new IllegalStateException("Empty matrix");
    }

    /** Return a ONE of the same Fraction type using x / x for a nonzero x. */
    private static Fraction anyOne(Fraction[][] H) {
        for (Fraction[] row : H) {
            for (Fraction f : row) {
                // find non-zero
                if (f.compareTo(f.subtract(f)) != 0) {
                    return f.divide(f);
                }
            }
        }
        // Fallback: 0/0 would be invalid; if all-zero rows exist that's degenerate.
        throw new IllegalStateException("Cannot construct 1 from all-zero matrix");
    }

    /** Check a0 + a·x >= 0 for all rows. */
    private static boolean feasible(Fraction[][] H, Fraction[] x, Fraction ZERO) {
        for (Fraction[] row : H) {
            Fraction lhs = row[0];
            for (int j = 0; j < x.length; j++) {
                lhs = lhs.add(row[j + 1].multiply(x[j]));
            }
            if (lhs.compareTo(ZERO) < 0) return false;
        }
        return true;
    }

    private static String canonical(Fraction[] v) {
        StringBuilder sb = new StringBuilder();
        for (Fraction f : v) sb.append(f.toString()).append('|');
        return sb.toString();
    }

    /** First k-combination [0..k-1]. */
    private static int[] initComb(int k) {
        int[] a = new int[k];
        for (int i = 0; i < k; i++) a[i] = i;
        return a;
    }

    /** Next k-combination of 0..n-1; null when done. */
    private static int[] nextComb(int[] a, int n, int k) {
        int i = k - 1;
        while (i >= 0 && a[i] == n - k + i) i--;
        if (i < 0) return null;
        a[i]++;
        for (int j = i + 1; j < k; j++) a[j] = a[j - 1] + 1;
        return a;
    }

    /**
     * Exact Gaussian elimination over Fractions for A x = b.
     * Uses add/subtract/multiply/divide; no negate()/zero()/one() needed.
     * Returns null if singular/inconsistent.
     */
    private static Fraction[] solve(Fraction[][] A, Fraction[] b, Fraction ZERO) {
        int n = b.length;
        Fraction[][] M = new Fraction[n][n + 1];
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, M[i], 0, n);
            M[i][n] = b[i];
        }

        int row = 0;
        for (int col = 0; col < n && row < n; col++) {
            int piv = row;
            while (piv < n && M[piv][col].compareTo(ZERO) == 0) piv++;
            if (piv == n) continue; // no pivot in this column

            if (piv != row) {
                Fraction[] tmp = M[piv]; M[piv] = M[row]; M[row] = tmp;
            }

            Fraction diag = M[row][col];
            if (diag.compareTo(ZERO) == 0) continue;

            // normalize pivot row
            for (int j = col; j <= n; j++) M[row][j] = M[row][j].divide(diag);

            // eliminate column from other rows
            for (int r = 0; r < n; r++) {
                if (r == row) continue;
                Fraction factor = M[r][col];
                if (factor.compareTo(ZERO) != 0) {
                    for (int j = col; j <= n; j++) {
                        M[r][j] = M[r][j].subtract(factor.multiply(M[row][j]));
                    }
                }
            }
            row++;
        }

        // inconsistency: [0 ... 0 | nonzero]
        for (int i = 0; i < n; i++) {
            boolean allZero = true;
            for (int j = 0; j < n; j++) {
                if (M[i][j].compareTo(ZERO) != 0) { allZero = false; break; }
            }
            if (allZero && M[i][n].compareTo(ZERO) != 0) return null;
        }

        // extract solution (assumes unique)
        Fraction[] x = new Fraction[n];
        for (int i = 0; i < n; i++) {
            int lead = -1;
            for (int j = 0; j < n; j++) {
                if (M[i][j].compareTo(ZERO) != 0) { lead = j; break; }
            }
            if (lead == -1) return null; // singular
            x[lead] = M[i][n];
        }
        for (int i = 0; i < n; i++) if (x[i] == null) return null;
        return x;
    }
}
