package com.vertexenumeration;

import java.util.*;

final class FacetEnumerator {

    private EnumStats lastStats = null;
    public EnumStats getLastStats() { return lastStats; }

    /** Main entry: convert V→H (supports rays). */
    public Polyhedron enumerate(Polyhedron vRep) {
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
            if (lead.compareTo(ZERO) == 0) rays.add(row);   // [0|r]
            else                           verts.add(row); // [1|x]
        }

        if (verts.size() + rays.size() < d) {
            this.lastStats = new EnumStats(); // no facets
            return new Polyhedron(Polyhedron.Type.H, 0, n, false, new Matrix(0, n));
        }

        List<Fraction[]> lifted = new ArrayList<>(verts.size() + rays.size());
        lifted.addAll(verts);
        lifted.addAll(rays);
        final int L = lifted.size();

        Map<String, Fraction[]> uniq = new HashMap<>();
        int[] comb = initComb(d);

        EnumStats stats = new EnumStats();

        while (comb != null) {
            Fraction[][] A = new Fraction[d][n];
            for (int i = 0; i < d; i++) {
                System.arraycopy(lifted.get(comb[i]), 0, A[i], 0, n);
            }

            Fraction[] h = nullspace1(A);
            if (h != null) {
                if (!isValidForAll(h, verts, rays)) {
                    Fraction[] neg = negate(h);
                    if (isValidForAll(neg, verts, rays)) h = neg; else h = null;
                }

                if (h != null) {
                    String key = canonicalFacetKey(h);
                    if (!uniq.containsKey(key)) {
                        uniq.put(key, canonicalizeRow(h));
                        stats.bases++;
                    }
                }
            }

            comb = nextComb(comb, L, d);
        }

        List<Fraction[]> facets = new ArrayList<>(uniq.values());

        stats.facets = facets.size();
        stats.bases  = facets.size();
        stats.setMode(LrsDat.Mode.CH);
        this.lastStats = stats;

        // ---- build matrix ----
        Matrix out = new Matrix(facets.size(), n);
        for (int i = 0; i < facets.size(); i++) {
            for (int j = 0; j < n; j++) out.set(i, j, facets.get(i)[j]);
        }

        return new Polyhedron(Polyhedron.Type.H, facets.size(), n, false, out);

    }

    // --- helpers ---
    private static boolean isValidForAll(Fraction[] h,
                                         List<Fraction[]> verts,
                                         List<Fraction[]> rays) {
        Fraction ZERO = h[0].subtract(h[0]);
        int n = h.length;
        for (Fraction[] row : verts) {
            Fraction s = ZERO;
            for (int j = 0; j < n; j++) s = s.add(row[j].multiply(h[j]));
            if (s.compareTo(ZERO) < 0) return false;
        }
        for (Fraction[] row : rays) {
            Fraction s = ZERO;
            for (int j = 0; j < n; j++) s = s.add(row[j].multiply(h[j]));
            if (s.compareTo(ZERO) < 0) return false;
        }
        return true;
    }

    private static Fraction[] canonicalizeRow(Fraction[] row) {
        int n = row.length, first = -1;
        Fraction ZERO = row[0].subtract(row[0]);
        for (int j = 0; j < n; j++) if (row[j].compareTo(ZERO) != 0) { first = j; break; }
        if (first == -1) return row.clone();
        Fraction s = row[first];
        if (s.compareTo(ZERO) < 0) s = ZERO.subtract(s);
        Fraction[] out = new Fraction[n];
        for (int j = 0; j < n; j++) out[j] = row[j].divide(s);
        return out;
    }

    private static String canonicalFacetKey(Fraction[] row) {
        return Arrays.toString(canonicalizeRow(row));
    }

    private static Fraction[] negate(Fraction[] v) {
        Fraction ZERO = v[0].subtract(v[0]);
        Fraction[] out = new Fraction[v.length];
        for (int i = 0; i < v.length; i++) out[i] = ZERO.subtract(v[i]);
        return out;
    }

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
                if (f.compareTo(ZERO) != 0)
                    for (int j = col; j < n; j++) M[i2][j] = M[i2][j].subtract(f.multiply(M[row][j]));
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
        Fraction[] x = new Fraction[n];
        Arrays.fill(x, ZERO);
        x[free] = Fraction.ONE;
        for (int i = 0; i < r; i++) if (lead[i] != -1) {
            int piv = lead[i];
            x[piv] = ZERO.subtract(M[i][free]);
        }
        return x;
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
}
