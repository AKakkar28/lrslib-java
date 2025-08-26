package com.vertexenumeration;

import java.util.*;

/** Exact facet enumerator for V-representations (polytopes).
 *  lrs-like ordering: facets sorted by lexicographically minimum cobasis.
 */
final class FacetEnumerator {

    private static final class Facet {
        final Fraction[] h;   // [a0, a1..ad]
        final int[] cobasis;  // size d, lex-min cobasis (vertex indices in input order)
        final String key;     // canonical dedup key
        Facet(Fraction[] h, int[] cobasis, String key){ this.h=h; this.cobasis=cobasis; this.key=key; }
    }

    static Polyhedron fromV(Polyhedron vRep) {
        final int m = vRep.getRowCount();
        final int n = vRep.getColCount(); // 1 + d
        final int d = n - 1;

        // ---- collect vertices (leading column == 1). Ignore rays for now. ----
        Matrix M = vRep.getMatrix();
        List<Fraction[]> verts = new ArrayList<>();
        for (int i = 0; i < m; i++) {
            Fraction lead = M.get(i, 0);
            Fraction ZERO = lead.subtract(lead);
            if (lead.compareTo(ZERO) == 0) continue;     // ray -> skip
            Fraction[] row = new Fraction[n];
            for (int j = 0; j < n; j++) row[j] = M.get(i, j);
            verts.add(row);
        }
        int mv = verts.size();
        if (mv < d + 1) {
            // Not enough points to define a d-polytope.
            return new Polyhedron(Polyhedron.Type.H, 0, n, false, new Matrix(0, n));
        }

        // ---- enumerate facets via all d-combinations of vertices ----
        Map<String,Facet> uniq = new HashMap<>();
        int[] comb = initComb(d);
        while (comb != null) {
            // Rows [1|x] of the chosen d vertices
            Fraction[][] A = new Fraction[d][n];
            for (int i = 0; i < d; i++) System.arraycopy(verts.get(comb[i]), 0, A[i], 0, n);

            // normal (a0,a) with A * (a0,a)^T = 0  (expect nullspace dim = 1)
            Fraction[] h = nullspace1(A);
            if (h != null) {
                // orient so all vertices satisfy >= 0
                if (!allGe(verts, h)) {
                    Fraction[] neg = negate(h);
                    if (allGe(verts, neg)) h = neg; else h = null;
                }
                if (h != null) {
                    String key = canonicalFacetKey(h);
                    if (!uniq.containsKey(key)) {
                        int[] cob = lexMinCobasis(h, verts, d);
                        uniq.put(key, new Facet(h, cob, key));
                    }
                }
            }

            comb = nextComb(comb, mv, d);
        }

        // ---- order facets by lrs rule: lexicographic cobasis ----
        List<Facet> facets = new ArrayList<>(uniq.values());
        facets.sort((p, q) -> lexInt(p.cobasis, q.cobasis));

        // ---- build output matrix in that order ----
        Matrix out = new Matrix(facets.size(), n);
        for (int i = 0; i < facets.size(); i++)
            for (int j = 0; j < n; j++)
                out.set(i, j, facets.get(i).h[j]);

        return new Polyhedron(Polyhedron.Type.H, facets.size(), n, false, out);
    }

    // ---------- lrs-style lex-min cobasis on tight set ----------

    // Return the lexicographically first size-d subset of tight vertices
    // (w.r.t. input order) that is affinely independent (rank d on rows [1|x]).
    private static int[] lexMinCobasis(Fraction[] h, List<Fraction[]> verts, int d) {
        // collect tight indices
        List<Integer> tight = new ArrayList<>();
        Fraction ZERO = h[0].subtract(h[0]);
        int n = h.length; // 1+d
        for (int i = 0; i < verts.size(); i++) {
            Fraction[] vi = verts.get(i); // [1, x...]
            Fraction s = ZERO;
            for (int j = 0; j < n; j++) s = s.add(vi[j].multiply(h[j]));
            if (s.compareTo(ZERO) == 0) tight.add(i);
        }
        // iterate d-combinations in lex order; first full-rank wins
        int t = tight.size();
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

        // fallback (shouldn’t happen for a valid facet)
        int[] cob = new int[d];
        for (int i = 0; i < d; i++) cob[i] = tight.get(i);
        return cob;
    }

    // rank of rows [1|x] by exact Gauss elimination
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
            for (int i = 0; i < r; i++) if (i != row) {
                Fraction f = M[i][col];
                if (f.compareTo(ZERO) != 0)
                    for (int j = col; j < n; j++) M[i][j] = M[i][j].subtract(f.multiply(M[row][j]));
            }
            row++;
        }
        return row;
    }

    // ---------- algebra helpers ----------

    private static Fraction[] negate(Fraction[] v){
        Fraction ZERO = v[0].subtract(v[0]);
        Fraction[] out = new Fraction[v.length];
        for (int i=0;i<v.length;i++) out[i] = ZERO.subtract(v[i]);
        return out;
    }

    /** Check a0 + a·x >= 0 for all vertices (rows are [1, x...]). */
    private static boolean allGe(List<Fraction[]> verts, Fraction[] h){
        Fraction ZERO = h[0].subtract(h[0]);
        int n = h.length;
        for (Fraction[] row : verts) {
            Fraction s = ZERO;
            for (int j = 0; j < n; j++) s = s.add(row[j].multiply(h[j]));
            if (s.compareTo(ZERO) < 0) return false;
        }
        return true;
    }

    /** Nullspace dimension 1 vector for A (r x n) where n = d+1, r = d. */
    private static Fraction[] nullspace1(Fraction[][] A){
        int r = A.length, n = A[0].length;
        Fraction ZERO = A[0][0].subtract(A[0][0]);
        Fraction[][] M = new Fraction[r][n];
        for (int i = 0; i < r; i++) System.arraycopy(A[i], 0, M[i], 0, n);

        int row = 0; int[] lead = new int[r]; Arrays.fill(lead, -1);
        for (int col = 0; col < n && row < r; col++){
            int p = row; while (p < r && M[p][col].compareTo(ZERO) == 0) p++;
            if (p == r) continue;
            if (p != row) { Fraction[] t = M[p]; M[p] = M[row]; M[row] = t; }
            Fraction diag = M[row][col];
            for (int j = col; j < n; j++) M[row][j] = M[row][j].divide(diag);
            for (int i = 0; i < r; i++) if (i != row){
                Fraction f = M[i][col];
                if (f.compareTo(ZERO) != 0)
                    for (int j = col; j < n; j++) M[i][j] = M[i][j].subtract(f.multiply(M[row][j]));
            }
            lead[row] = col; row++;
        }

        int rank = 0; for (int i=0;i<r;i++) if (lead[i]!=-1) rank++;
        if (n - rank != 1) return null;

        boolean[] isPivot = new boolean[n];
        for (int i=0;i<r;i++) if (lead[i]>=0) isPivot[lead[i]] = true;

        int free = -1; for (int j = n-1; j>=0; j--) if (!isPivot[j]) { free=j; break; }
        if (free == -1) return null;

        Fraction ONE = null;
        outer:
        for (int i=0;i<r;i++) for (int j=0;j<n;j++)
            if (M[i][j].compareTo(ZERO)!=0){ ONE = M[i][j].divide(M[i][j]); break outer; }
        if (ONE == null) return null;

        Fraction[] x = new Fraction[n];
        for (int j=0;j<n;j++) x[j]=ZERO;
        x[free]=ONE;
        for (int i=0;i<r;i++) if (lead[i]!=-1){
            int piv = lead[i];
            x[piv] = ZERO.subtract(M[i][free]);
        }
        return x;
    }

    /** Canonicalize [a0,a...] so first nonzero is +1; use as dedup key. */
    private static String canonicalFacetKey(Fraction[] row){
        int n=row.length, first=-1;
        Fraction ZERO = row[0].subtract(row[0]);
        for (int j=0;j<n;j++) if (row[j].compareTo(ZERO)!=0){ first=j; break; }
        if (first==-1) return "0";
        Fraction s = row[first];
        if (s.compareTo(ZERO) < 0) s = ZERO.subtract(s);
        StringBuilder sb = new StringBuilder();
        for (int j=0;j<n;j++) sb.append(row[j].divide(s).toString()).append('|');
        return sb.toString();
    }

    // ---------- small utilities ----------

    private static int[] initComb(int k){ int[] a=new int[k]; for(int i=0;i<k;i++) a[i]=i; return a; }
    private static int[] nextComb(int[] a, int n, int k){
        int i = k - 1;
        while (i >= 0 && a[i] == n - k + i) i--;
        if (i < 0) return null;
        a[i]++;
        for (int j = i + 1; j < k; j++) a[j] = a[j - 1] + 1;
        return a;
    }
    private static boolean nextCombInPlace(int[] a, int n, int k){
        int i = k - 1;
        while (i >= 0 && a[i] == n - k + i) i--;
        if (i < 0) return false;
        a[i]++;
        for (int j = i + 1; j < k; j++) a[j] = a[j - 1] + 1;
        return true;
    }
    private static int lexInt(int[] A, int[] B){
        for (int i=0;i<A.length;i++){ int c = Integer.compare(A[i], B[i]); if (c!=0) return c; }
        return 0;
    }
}
