package com.vertexenumeration;

import java.util.*;

/** Exact facet enumerator for V-representations (vertices + rays).
 *  - Accepts rows [1 | x] (vertices) and [0 | r] (rays).
 *  - Enumerates facets by testing all d-combinations of lifted rows.
 *  - Valid facet h=(a0,a):  a0 + a·x ≥ 0 for all vertices and a·r ≥ 0 for all rays.
 *  - Deduplicates by canonical scaling (first nonzero = +1).
 *  - Ordering (temporary but close to lrs):
 *        (1) a0 == 0 first,
 *        (2) lex-min cobasis of tight *vertices* (input order),
 *        (3) canonical row string as tiebreak.
 */
final class FacetEnumerator {

    /** Internal record to sort and print facets. */
    private static final class Facet {
        final Fraction[] h;     // [a0, a1..ad]
        final int[] cobasis;    // size ≤ d; lex-min tight vertex indices (input order)
        final boolean a0isZero;
        final String tie;       // canonical row string for final tie-break
        final String key;       // canonical dedup key (first nonzero scaled to +1)
        Facet(Fraction[] h, int[] cobasis, boolean a0isZero, String tie, String key) {
            this.h = h; this.cobasis = cobasis; this.a0isZero = a0isZero; this.tie = tie; this.key = key;
        }
    }

    // Add a uniform entrypoint for V->H
    public Polyhedron enumerate(Polyhedron V) {
        Polyhedron out = fromV(V);   // <-- delegate to your existing method name
        return out;
    }

    private EnumStats lastStats = null;
    public EnumStats getLastStats() { return lastStats; }


    /** Convert V→H (supports rays). */
    static Polyhedron fromV(Polyhedron vRep) {
        final int m = vRep.getRowCount();
        final int n = vRep.getColCount();  // n = 1 + d
        final int d = n - 1;

        Matrix M = vRep.getMatrix();

        // Split input rows into vertices [1|x] and rays [0|r]
        List<Fraction[]> verts = new ArrayList<>();
        List<Fraction[]> rays  = new ArrayList<>();
        for (int i = 0; i < m; i++) {
            Fraction lead = M.get(i, 0);
            Fraction ZERO = lead.subtract(lead);
            Fraction[] row = new Fraction[n];
            for (int j = 0; j < n; j++) row[j] = M.get(i, j);
            if (lead.compareTo(ZERO) == 0) rays.add(row);  // [0|r]
            else                             verts.add(row); // [1|x]
        }

        // If total lifted rows < d, we cannot define a facet
        if (verts.size() + rays.size() < d) {
            return new Polyhedron(Polyhedron.Type.H, 0, n, false, new Matrix(0, n));
        }

        // Build unified lifted list (choose d rows from it)
        List<Fraction[]> lifted = new ArrayList<>(verts.size() + rays.size());
        lifted.addAll(verts);
        lifted.addAll(rays);
        final int L = lifted.size();

        // Enumerate all d-combinations to get candidate facets
        Map<String, Facet> uniq = new HashMap<>();
        int[] comb = initComb(d);
        while (comb != null) {
            // Build a d x (d+1) matrix A with chosen lifted rows
            Fraction[][] A = new Fraction[d][n];
            for (int i = 0; i < d; i++) {
                System.arraycopy(lifted.get(comb[i]), 0, A[i], 0, n);
            }

            // Find nullspace vector h = (a0, a) with A * h = 0 (expect dim 1)
            Fraction[] h = nullspace1(A);
            if (h != null) {
                // Orient so that all vertices/rays satisfy the facet inequality
                if (!isValidForAll(h, verts, rays)) {
                    Fraction[] neg = negate(h);
                    if (isValidForAll(neg, verts, rays)) h = neg; else h = null;
                }

                if (h != null) {
                    // Dedup canonically (first nonzero scaled to +1)
                    String key = canonicalFacetKey(h);
                    if (!uniq.containsKey(key)) {
                        // lrs reports cobasis using *vertices* only; pick lex-min tight vertex set of size ≤ d
                        int[] cob = lexMinCobasisVerticesOnly(h, verts, d);
                        boolean a0isZero = isZero(h[0]);
                        String tie = canonicalRow(h);
                        uniq.put(key, new Facet(h, cob, a0isZero, tie, key));
                    }
                }
            }

            comb = nextComb(comb, L, d);
        }

        // Sort facets: (a0 == 0 first) → lex-min cobasis → canonical row
        List<Facet> facets = new ArrayList<>(uniq.values());
        facets.sort((p, q) -> {
            if (p.a0isZero != q.a0isZero) return p.a0isZero ? -1 : 1;
            int c = lexInt(p.cobasis, q.cobasis);
            if (c != 0) return c;
            return p.tie.compareTo(q.tie);
        });

        // Build H matrix in that order
        Matrix out = new Matrix(facets.size(), n);
        for (int i = 0; i < facets.size(); i++) {
            for (int j = 0; j < n; j++) out.set(i, j, facets.get(i).h[j]);
        }
        return new Polyhedron(Polyhedron.Type.H, facets.size(), n, false, out);
    }

    // ---------- validity & cobasis ----------

    /** Check a0 + a·x ≥ 0 for all vertices and a·r ≥ 0 for all rays. */
    private static boolean isValidForAll(Fraction[] h,
                                         List<Fraction[]> verts,
                                         List<Fraction[]> rays) {
        Fraction ZERO = h[0].subtract(h[0]);
        int n = h.length;

        // vertices [1|x]
        for (Fraction[] row : verts) {
            Fraction s = ZERO;
            for (int j = 0; j < n; j++) s = s.add(row[j].multiply(h[j]));
            if (s.compareTo(ZERO) < 0) return false;
        }
        // rays [0|r] → same dot = a·r
        for (Fraction[] row : rays) {
            Fraction s = ZERO;
            for (int j = 0; j < n; j++) s = s.add(row[j].multiply(h[j]));
            if (s.compareTo(ZERO) < 0) return false;
        }
        return true;
    }

    /** Lex-first size-d subset of tight *vertices* (input order) with rank d on [1|x].
     *  If fewer than d tight vertices exist (unbounded facet), return as many as available. */
    private static int[] lexMinCobasisVerticesOnly(Fraction[] h, List<Fraction[]> verts, int d) {
        List<Integer> tight = new ArrayList<>();
        Fraction ZERO = h[0].subtract(h[0]);
        int n = h.length; // 1 + d
        for (int i = 0; i < verts.size(); i++) {
            Fraction[] vi = verts.get(i); // [1, x...]
            Fraction s = ZERO;
            for (int j = 0; j < n; j++) s = s.add(vi[j].multiply(h[j]));
            if (s.compareTo(ZERO) == 0) tight.add(i);
        }

        int t = tight.size();
        if (t < d) {
            // Not enough tight vertices: facet supported also by rays.
            // Return the lex-first subset that exists (size t).
            int size = Math.max(0, Math.min(d, t));
            int[] cobSmall = new int[size];
            for (int i = 0; i < size; i++) cobSmall[i] = tight.get(i);
            return cobSmall;
        }

        // Try all d-combinations of tight vertices (lex order) until rank d
        int[] a = new int[d];
        for (int i = 0; i < d; i++) a[i] = i;
        do {
            List<Fraction[]> rows = new ArrayList<>(d);
            for (int k = 0; k < d; k++) rows.add(verts.get(tight.get(a[k])));
            if (rank(rows) == d) {
                int[] cob = new int[d];
                for (int k = 0; k < d; k++) cob[k] = tight.get(a[k]);
                return cob;
            }
        } while (nextCombInPlace(a, t, d));

        // Fallback: first d tight vertices
        int[] cob = new int[d];
        for (int i = 0; i < d; i++) cob[i] = tight.get(i);
        return cob;
    }

    /** Rank of rows [1|x] with exact Gauss elimination. */
    private static int rank(List<Fraction[]> rows) {
        if (rows.isEmpty()) return 0;
        int r = rows.size(), n = rows.get(0).length;
        Fraction[][] M = new Fraction[r][n];
        for (int i = 0; i < r; i++) System.arraycopy(rows.get(i), 0, M[i], 0, n);
        Fraction ZERO = M[0][0].subtract(M[0][0]);

        int row = 0;
        for (int col = 0; col < n && row < r; col++) {
            int p = row;
            while (p < r && M[p][col].compareTo(ZERO) == 0) p++;
            if (p == r) continue;
            if (p != row) { Fraction[] t = M[p]; M[p] = M[row]; M[row] = t; }
            Fraction diag = M[row][col];
            for (int j = col; j < n; j++) M[row][j] = M[row][j].divide(diag);
            for (int i2 = 0; i2 < r; i2++) if (i2 != row) {
                Fraction f = M[i2][col];
                if (f.compareTo(ZERO) != 0)
                    for (int j = col; j < n; j++) M[i2][j] = M[i2][j].subtract(f.multiply(M[row][j]));
            }
            row++;
        }
        return row;
    }

    // ---------- algebra helpers ----------

    private static boolean isZero(Fraction f) { return f.compareTo(f.subtract(f)) == 0; }

    private static Fraction[] negate(Fraction[] v) {
        Fraction ZERO = v[0].subtract(v[0]);
        Fraction[] out = new Fraction[v.length];
        for (int i = 0; i < v.length; i++) out[i] = ZERO.subtract(v[i]);
        return out;
    }

    /** Nullspace dimension-1 vector for A (r × n), with n = d+1 and r = d. */
    private static Fraction[] nullspace1(Fraction[][] A) {
        int r = A.length, n = A[0].length;
        Fraction ZERO = A[0][0].subtract(A[0][0]);

        Fraction[][] M = new Fraction[r][n];
        for (int i = 0; i < r; i++) System.arraycopy(A[i], 0, M[i], 0, n);

        int row = 0;
        int[] lead = new int[r];
        Arrays.fill(lead, -1);

        for (int col = 0; col < n && row < r; col++) {
            int p = row;
            while (p < r && M[p][col].compareTo(ZERO) == 0) p++;
            if (p == r) continue;
            if (p != row) { Fraction[] t = M[p]; M[p] = M[row]; M[row] = t; }
            Fraction diag = M[row][col];
            for (int j = col; j < n; j++) M[row][j] = M[row][j].divide(diag);
            for (int i2 = 0; i2 < r; i2++) if (i2 != row) {
                Fraction f = M[i2][col];
                if (f.compareTo(ZERO) != 0) {
                    for (int j = col; j < n; j++) M[i2][j] = M[i2][j].subtract(f.multiply(M[row][j]));
                }
            }
            lead[row] = col;
            row++;
        }

        int rank = 0;
        for (int i = 0; i < r; i++) if (lead[i] != -1) rank++;
        if (n - rank != 1) return null;

        boolean[] isPivot = new boolean[n];
        for (int i = 0; i < r; i++) if (lead[i] >= 0) isPivot[lead[i]] = true;

        int free = -1;
        for (int j = n - 1; j >= 0; j--) if (!isPivot[j]) { free = j; break; }
        if (free == -1) return null;

        // Build ONE safely from any nonzero
        Fraction ONE = null;
        outer:
        for (int i = 0; i < r; i++)
            for (int j = 0; j < n; j++)
                if (M[i][j].compareTo(ZERO) != 0) { ONE = M[i][j].divide(M[i][j]); break outer; }
        if (ONE == null) return null;

        Fraction[] x = new Fraction[n];
        for (int j = 0; j < n; j++) x[j] = ZERO;
        x[free] = ONE;
        for (int i = 0; i < r; i++) if (lead[i] != -1) {
            int piv = lead[i];
            x[piv] = ZERO.subtract(M[i][free]);
        }
        return x;
    }

    /** Canonicalize [a0, a...] so the first nonzero equals +1; stringify for dedup. */
    private static String canonicalFacetKey(Fraction[] row) {
        int n = row.length, first = -1;
        Fraction ZERO = row[0].subtract(row[0]);
        for (int j = 0; j < n; j++) if (row[j].compareTo(ZERO) != 0) { first = j; break; }
        if (first == -1) return "0";
        Fraction s = row[first];
        if (s.compareTo(ZERO) < 0) s = ZERO.subtract(s);
        StringBuilder sb = new StringBuilder(n * 8);
        for (int j = 0; j < n; j++) sb.append(row[j].divide(s).toString()).append('|');
        return sb.toString();
    }

    /** Raw canonical string (no sign flip); used only as a sort tiebreak. */
    private static String canonicalRow(Fraction[] row) {
        StringBuilder sb = new StringBuilder(row.length * 8);
        for (Fraction f : row) sb.append(f.toString()).append('|');
        return sb.toString();
    }

    // ---------- small combinations & compares ----------

    private static int[] initComb(int k) { int[] a = new int[k]; for (int i = 0; i < k; i++) a[i] = i; return a; }

    private static int[] nextComb(int[] a, int n, int k) {
        int i = k - 1;
        while (i >= 0 && a[i] == n - k + i) i--;
        if (i < 0) return null;
        a[i]++;
        for (int j = i + 1; j < k; j++) a[j] = a[j - 1] + 1;
        return a;
    }

    private static boolean nextCombInPlace(int[] a, int n, int k) {
        int i = k - 1;
        while (i >= 0 && a[i] == n - k + i) i--;
        if (i < 0) return false;
        a[i]++;
        for (int j = i + 1; j < k; j++) a[j] = a[j - 1] + 1;
        return true;
    }

    private static int lexInt(int[] A, int[] B) {
        int la = (A == null) ? 0 : A.length;
        int lb = (B == null) ? 0 : B.length;
        int L = Math.min(la, lb);
        for (int i = 0; i < L; i++) {
            int c = Integer.compare(A[i], B[i]);
            if (c != 0) return c;
        }
        // shorter list sorts first if equal prefix
        return Integer.compare(la, lb);
    }
}
