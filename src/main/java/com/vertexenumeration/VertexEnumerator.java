package com.vertexenumeration;

import java.util.*;

public class VertexEnumerator {

    private EnumStats lastStats = null;
    public EnumStats getLastStats() { return lastStats; }

    public Polyhedron enumerate(Polyhedron input) {
        if (input.getType() == Polyhedron.Type.H) {
            return enumerateFromH(input);   // H → V
        } else {
            // V → H
            FacetEnumerator enumerator = new FacetEnumerator();
            Polyhedron out = enumerator.enumerate(input);

            lastStats = enumerator.getLastStats();
            lastStats.vertices = 0;
            lastStats.rays = 0;
            lastStats.bases = out.getRowCount();
            lastStats.integerVertices = 0;
            lastStats.maxDepth = 0;

            return out;
        }
    }

    private Polyhedron enumerateFromH(Polyhedron hRep) {
        final int m = hRep.getRowCount();
        final int n = hRep.getColCount(); // 1 + dim
        final int d = n - 1;

        if (d <= 0 || m < d) {
            return new Polyhedron(Polyhedron.Type.V, 0, n, false, new Matrix(0, n));
        }

        // Convert H-polyhedron into matrix
        Matrix Mm = hRep.getMatrix();
        Fraction[][] H = new Fraction[m][n];
        for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) H[i][j] = Mm.get(i, j);

        // ---- delegate to ReverseSearchEnumerator ----
        ReverseSearchEnumerator.Result res = new ReverseSearchEnumerator(H).run();

        this.lastStats = res.stats;

        // ---- rays ----
        List<RayRec> rays = enumerateExtremeRaysWithTightSets(H);
        lastStats.rays = rays.size();

        // ---- vertices from reverse search ----
        List<Fraction[]> verts = res.vertices;
        List<int[]> vertexBases = new ArrayList<>(); // TODO: carry bases if needed

        // ---- attach rays ----
        Map<Integer, List<Fraction[]>> raysAt = new HashMap<>();
        for (int i = 0; i < verts.size(); i++) raysAt.put(i, new ArrayList<>());
        for (RayRec rr : rays) {
            int at = attachIndexForRay(H, vertexBases, rr.ray);
            if (at < 0) at = 0;
            raysAt.get(at).add(rr.ray);
        }

        // ---- build V matrix ----
        int totalRows = verts.size();
        for (List<Fraction[]> lst : raysAt.values()) totalRows += lst.size();

        Matrix V = new Matrix(totalRows, n);
        int row = 0;
        for (int i = 0; i < verts.size(); i++) {
            Fraction[] v = verts.get(i);
            for (int j = 0; j < n; j++) V.set(row, j, v[j]);
            row++;
            List<Fraction[]> here = raysAt.get(i);
            if (here != null) {
                for (Fraction[] r : here) {
                    for (int j = 0; j < n; j++) V.set(row, j, r[j]);
                    row++;
                }
            }
        }
        lastStats.setMode(LrsDat.Mode.VE);

        return new Polyhedron(Polyhedron.Type.V, totalRows, n, false, V);
    }

    // ---------- Support classes / methods ----------

    static final class RayRec {
        final Fraction[] ray;
        final int[] tight;
        final String key;
        RayRec(Fraction[] ray, int[] tight, String key) {
            this.ray = ray; this.tight = tight; this.key = key;
        }
    }

    private static boolean isIntegerVertex(Fraction[] v) {
        for (int j = 1; j < v.length; j++) {
            if (!v[j].denominator().equals(java.math.BigInteger.ONE)) return false;
        }
        return true;
    }

    private static int attachIndexForRay(Fraction[][] H, List<int[]> vertexBases, Fraction[] ray) {
        final int d = H[0].length - 1;
        Fraction[] v = new Fraction[d];
        System.arraycopy(ray, 1, v, 0, d);

        for (int i = 0; i < vertexBases.size(); i++) {
            int zeros = 0;
            Fraction nonZero = null;
            for (int rIdx : vertexBases.get(i)) {
                Fraction s = dotRowA(H, rIdx, v);
                if (s.compareTo(Fraction.ZERO) == 0) {
                    zeros++;
                } else {
                    nonZero = s;
                }
            }
            if (zeros == d - 1 && (nonZero == null || nonZero.compareTo(Fraction.ZERO) > 0)) {
                return i;
            }
        }
        return -1;
    }

    private static Fraction dotRowA(Fraction[][] H, int row, Fraction[] v) {
        int d = v.length;
        Fraction s = Fraction.ZERO;
        for (int j = 0; j < d; j++) s = s.add(H[row][j + 1].multiply(v[j]));
        return s;
    }

    private List<RayRec> enumerateExtremeRaysWithTightSets(Fraction[][] H) {
        final int m = H.length, n = H[0].length, d = n - 1;
        if (d - 1 <= 0) return Collections.emptyList();
        Fraction[][] A = new Fraction[m][d];
        for (int i = 0; i < m; i++) for (int j = 0; j < d; j++) A[i][j] = H[i][j+1];

        List<RayRec> out = new ArrayList<>();
        Set<String> seen = new HashSet<>();

        int[] comb = initComb(d - 1);
        while (comb != null) {
            Fraction[][] As = new Fraction[d - 1][d];
            for (int i = 0; i < d - 1; i++) {
                int r = comb[i];
                for (int j = 0; j < d; j++) As[i][j] = A[r][j];
            }

            Fraction[] v = nullspace1(As);
            if (v != null) {
                if (!allGe(A, v)) {
                    Fraction[] vneg = negate(v);
                    if (allGe(A, vneg)) v = vneg; else v = null;
                }
                if (v != null) {
                    Fraction[] ray = new Fraction[d + 1];
                    ray[0] = Fraction.ZERO;
                    System.arraycopy(v, 0, ray, 1, d);
                    String rkey = rayKey(ray);
                    if (seen.add(rkey)) out.add(new RayRec(ray, comb.clone(), rkey));
                }
            }
            comb = nextComb(comb, m, d - 1);
        }
        return out;
    }

    private static boolean allGe(Fraction[][] A, Fraction[] v) {
        for (int i = 0; i < A.length; i++) {
            Fraction s = Fraction.ZERO;
            for (int j = 0; j < v.length; j++) s = s.add(A[i][j].multiply(v[j]));
            if (s.compareTo(Fraction.ZERO) < 0) return false;
        }
        return true;
    }

    private static Fraction[] nullspace1(Fraction[][] As) {
        int r = As.length, d = As[0].length;
        Fraction[][] M = new Fraction[r][d];
        for (int i = 0; i < r; i++) System.arraycopy(As[i], 0, M[i], 0, d);
        int row = 0;
        int[] lead = new int[r];
        Arrays.fill(lead, -1);
        for (int col = 0; col < d && row < r; col++) {
            int p = row;
            while (p < r && M[p][col].compareTo(Fraction.ZERO) == 0) p++;
            if (p == r) continue;
            if (p != row) { Fraction[] t = M[p]; M[p] = M[row]; M[row] = t; }
            Fraction diag = M[row][col];
            for (int j = col; j < d; j++) M[row][j] = M[row][j].divide(diag);
            for (int i = 0; i < r; i++) if (i != row) {
                Fraction f = M[i][col];
                if (f.compareTo(Fraction.ZERO) != 0) {
                    for (int j = col; j < d; j++) {
                        M[i][j] = M[i][j].subtract(f.multiply(M[row][j]));
                    }
                }
            }
            lead[row] = col;
            row++;
        }
        int rank = 0;
        for (int i = 0; i < r; i++) if (lead[i] != -1) rank++;
        if (d - rank != 1) return null;
        boolean[] isPivotCol = new boolean[d];
        for (int i = 0; i < r; i++) if (lead[i] >= 0) isPivotCol[lead[i]] = true;
        int free = -1;
        for (int j = d - 1; j >= 0; j--) if (!isPivotCol[j]) { free = j; break; }
        if (free == -1) return null;
        Fraction[] v = new Fraction[d];
        Arrays.fill(v, Fraction.ZERO);
        v[free] = Fraction.ONE;
        for (int i = 0; i < r; i++) if (lead[i] != -1) {
            int piv = lead[i];
            v[piv] = Fraction.ZERO.subtract(M[i][free]);
        }
        return v;
    }

    private static Fraction[] negate(Fraction[] v){
        Fraction[] o = new Fraction[v.length];
        for (int i=0;i<v.length;i++) o[i] = Fraction.ZERO.subtract(v[i]);
        return o;
    }

    private static String rayKey(Fraction[] ray) {
        int n = ray.length;
        int first = -1;
        for (int j = 1; j < n; j++) if (ray[j].compareTo(Fraction.ZERO) != 0) { first = j; break; }
        if (first == -1) first = 1;
        Fraction s = ray[first];
        if (s.compareTo(Fraction.ZERO) < 0) s = Fraction.ZERO.subtract(s);
        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < n; j++) {
            Fraction rj = (j == 0) ? ray[0] : ray[j].divide(s);
            sb.append(rj.toString()).append('|');
        }
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
}
