package com.vertexenumeration;

import java.util.*;

/** Minimal exact facet enumerator for V-representations (polytopes only).
 *  Given m rows [1 v1 ... vd], enumerates supporting hyperplanes a0 + a·x >= 0
 *  by testing all d-combinations of vertices and keeping valid, deduped facets.
 *  Rays in the input are ignored for now (TODO).
 */
final class FacetEnumerator {

    static Polyhedron fromV(Polyhedron vRep) {
        final int m = vRep.getRowCount();
        final int n = vRep.getColCount();      // 1 + d
        final int d = n - 1;

        Matrix M = vRep.getMatrix();

        // Collect only vertices (leading 1). Ignore rays (leading 0) for now.
        List<Fraction[]> verts = new ArrayList<>();
        for (int i = 0; i < m; i++) {
            Fraction lead = M.get(i, 0);
            Fraction ZERO = lead.subtract(lead);
            if (lead.compareTo(ZERO) == 0) continue; // a ray; skip in this minimal pass
            Fraction[] row = new Fraction[n];
            for (int j = 0; j < n; j++) row[j] = M.get(i, j);
            verts.add(row);
        }
        int mv = verts.size();
        if (mv < d+1) {
            // Degenerate or not enough points; return empty H for now.
            return new Polyhedron(Polyhedron.Type.H, 0, n, false, new Matrix(0, n));
        }

        // Enumerate all d-combinations of vertices to get a candidate facet.
        Set<String> seen = new HashSet<>();
        List<Fraction[]> facets = new ArrayList<>();

        int[] comb = initComb(d);
        while (comb != null) {
            // Build A: rows are [1 | v_i], size d x (d+1); find (a0,a) s.t. A*(a0,a) = 0
            Fraction[][] A = new Fraction[d][n];
            for (int i = 0; i < d; i++) {
                Fraction[] vi = verts.get(comb[i]);
                // vi already contains [1, v...]; copy directly
                System.arraycopy(vi, 0, A[i], 0, n);
            }

            Fraction[] normal = nullspace1(A);   // returns (a0,a) or null if nullspace dim != 1
            if (normal != null) {
                // Orient so that all vertices satisfy a0 + a·x >= 0
                if (!allGe(verts, normal)) {
                    Fraction[] neg = negate(normal);
                    if (allGe(verts, neg)) normal = neg;
                    else normal = null; // not a supporting hyperplane
                }

                if (normal != null) {
                    String key = canonicalFacetKey(normal);
                    if (seen.add(key)) {
                        // Store as H-row: [a0, a...]
                        facets.add(normal);
                    }
                }
            }

            comb = nextComb(comb, mv, d);
        }

        // Build output matrix
        Matrix H = new Matrix(facets.size(), n);
        for (int i = 0; i < facets.size(); i++) {
            for (int j = 0; j < n; j++) H.set(i, j, facets.get(i)[j]);
        }
        return new Polyhedron(Polyhedron.Type.H, facets.size(), n, false, H);
    }

    // ---------- helpers ----------

    private static int[] initComb(int k){ int[] a=new int[k]; for(int i=0;i<k;i++) a[i]=i; return a; }
    private static int[] nextComb(int[] a, int n, int k){
        int i = k - 1;
        while (i >= 0 && a[i] == n - k + i) i--;
        if (i < 0) return null;
        a[i]++;
        for (int j = i + 1; j < k; j++) a[j] = a[j - 1] + 1;
        return a;
    }

    private static Fraction[] negate(Fraction[] v){
        Fraction ZERO = v[0].subtract(v[0]);
        Fraction[] out = new Fraction[v.length];
        for (int i=0;i<v.length;i++) out[i] = ZERO.subtract(v[i]);
        return out;
    }

    /** Check a0 + a·x >= 0 for all vertices (rows are [1, x...]). */
    private static boolean allGe(List<Fraction[]> verts, Fraction[] normal){
        Fraction ZERO = normal[0].subtract(normal[0]);
        int d = normal.length - 1;
        for (Fraction[] row : verts) {
            Fraction s = ZERO;
            for (int j = 0; j <= d; j++) s = s.add(row[j].multiply(normal[j])); // row = [1, x...]
            if (s.compareTo(ZERO) < 0) return false;
        }
        return true;
    }

    /** Find one nonzero vector in nullspace( A ) for A (r x n), expecting n = r+1 so nullspace is 1D. */
    private static Fraction[] nullspace1(Fraction[][] A){
        int r = A.length, n = A[0].length;
        Fraction ZERO = A[0][0].subtract(A[0][0]);

        Fraction[][] M = new Fraction[r][n];
        for (int i = 0; i < r; i++) System.arraycopy(A[i], 0, M[i], 0, n);

        int row = 0;
        int[] lead = new int[r];
        Arrays.fill(lead, -1);

        for (int col = 0; col < n && row < r; col++){
            int p = row;
            while (p < r && M[p][col].compareTo(ZERO) == 0) p++;
            if (p == r) continue;
            if (p != row) { Fraction[] t = M[p]; M[p] = M[row]; M[row] = t; }
            Fraction diag = M[row][col];
            for (int j = col; j < n; j++) M[row][j] = M[row][j].divide(diag);
            for (int i = 0; i < r; i++) if (i != row){
                Fraction f = M[i][col];
                if (f.compareTo(ZERO) != 0){
                    for (int j = col; j < n; j++) M[i][j] = M[i][j].subtract(f.multiply(M[row][j]));
                }
            }
            lead[row] = col;
            row++;
        }

        int rank = 0;
        for (int i = 0; i < r; i++) if (lead[i] != -1) rank++;
        if (n - rank != 1) return null;

        // pick a free column and set it to 1, then back-substitute to get a null vector
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
        if (ONE == null) return null; // A==0 -> not a facet in this context

        Fraction[] x = new Fraction[n];
        for (int j = 0; j < n; j++) x[j] = ZERO;
        x[free] = ONE;
        for (int i = 0; i < r; i++) if (lead[i] != -1) {
            int piv = lead[i];
            x[piv] = ZERO.subtract(M[i][free]);
        }
        return x;
    }

    /** Canonicalize facet row [a0, a...] so first nonzero is +1; use as a dedup key. */
    private static String canonicalFacetKey(Fraction[] row){
        int n = row.length;
        int first = -1;
        Fraction ZERO = row[0].subtract(row[0]);
        for (int j = 0; j < n; j++) if (row[j].compareTo(ZERO) != 0) { first = j; break; }
        if (first == -1) return "0"; // all zero shouldn't happen here

        Fraction s = row[first];
        if (s.compareTo(ZERO) < 0) s = ZERO.subtract(s); // make it positive
        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < n; j++) sb.append(row[j].divide(s).toString()).append('|');
        return sb.toString();
    }
}
