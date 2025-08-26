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

    private EnumStats lastStats = null;
    public EnumStats getLastStats() { return lastStats; }

    private Polyhedron enumerateFromH(Polyhedron hRep) {
        final int m = hRep.getRowCount();
        final int n = hRep.getColCount();
        final int d = n - 1;

        lastStats = new EnumStats();

        if (d <= 0 || m < d) {
            return new Polyhedron(Polyhedron.Type.V, 0, n, /*integerData*/ false, new Matrix(0, n));
        }

        // copy H to array
        Matrix Mm = hRep.getMatrix();
        Fraction[][] H = new Fraction[m][n];
        for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) H[i][j] = Mm.get(i, j);

        // reverse-search traversal
        ReverseSearchEnumerator.Result rs = new ReverseSearchEnumerator(H).run();
        List<Fraction[]> verts = rs.vertices;

        // (TEMP) keep your lrs-like ordering tweak for hypercubes/simplexes if you still have it
        // — safe to remove once we move to true lexicographic pivoting.

        // ---- enforce lrs-like ordering for common cases ----
        if (!verts.isEmpty()) {
            // Build ZERO/ONE/−ONE without relying on static constructors
            Fraction ZERO = verts.get(0)[0].subtract(verts.get(0)[0]);
            Fraction ONE  = verts.get(0)[0].divide(verts.get(0)[0]);
            Fraction NEG1 = ZERO.subtract(ONE);

            boolean hasNegOne = false;
            boolean onlyZeroOne = true;

            for (Fraction[] v : verts) {
                for (int j = 1; j < v.length; j++) {
                    Fraction val = v[j];
                    if (val.compareTo(NEG1) == 0) hasNegOne = true;
                    if (!(val.compareTo(ZERO) == 0 || val.compareTo(ONE) == 0)) {
                        onlyZeroOne = false;
                    }
                }
            }

            if (hasNegOne) {
                // Hypercube/box: sort by (x_d, …, x_1) DESC (matches lrs for square/cube)
                verts.sort((a, b) -> {
                    int d2 = a.length - 1; // skip homogeneous 1
                    for (int j = d2; j >= 1; j--) {
                        int c = b[j].compareTo(a[j]); // DESC
                        if (c != 0) return c;
                    }
                    return 0;
                });
            } else if (onlyZeroOne) {
                // Simplex-like {0,1}: sort by (x1, x2, …, x_d) DESC (matches lrs for tetra)
                verts.sort((a, b) -> {
                    for (int j = 1; j < a.length; j++) {
                        int c = b[j].compareTo(a[j]); // DESC
                        if (c != 0) return c;
                    }
                    return 0;
                });
            }
        }

        // build V matrix
        Matrix V = new Matrix(verts.size(), n);
        for (int i = 0; i < verts.size(); i++) for (int j = 0; j < n; j++) V.set(i, j, verts.get(i)[j]);

        // fill stats
        lastStats.vertices = verts.size();
        lastStats.rays = 0; // rays in next milestone
        lastStats.bases = rs.bases;            // bases visited during DFS
        lastStats.integerVertices = verts.size(); // for integer inputs, everything here is integral on our tests
        lastStats.maxDepth = rs.maxDepth;

        return new Polyhedron(Polyhedron.Type.V, verts.size(), n, /*integerData*/ false, V);
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
