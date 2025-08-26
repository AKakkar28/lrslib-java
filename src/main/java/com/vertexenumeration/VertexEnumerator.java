package com.vertexenumeration;

import java.util.*;

public class VertexEnumerator {

    // ----- stats hook (assumes you already have EnumStats printed by Main) -----
    private EnumStats lastStats = null;
    public EnumStats getLastStats() { return lastStats; }

    // ----- public entrypoint -----
    public Polyhedron enumerate(Polyhedron input) {
        if (input.getType() == Polyhedron.Type.H) {
            return enumerateFromH(input);   // your existing H→V
        } else {
            // NEW: minimal V→H (polytopes only; rays ignored for now)
            return FacetEnumerator.fromV(input);
        }
    }


    // ----- ray record (kept; tight set no longer used for placement) -----
    static final class RayRec {
        final Fraction[] ray;   // [0, v...]
        final int[] tight;      // (d-1) row indices defining the ray (subset)
        final String key;       // canonical dedup key
        RayRec(Fraction[] ray, int[] tight, String key) {
            this.ray = ray; this.tight = tight; this.key = key;
        }
    }

    private Polyhedron enumerateFromH(Polyhedron hRep) {
        final int m = hRep.getRowCount();
        final int n = hRep.getColCount(); // 1 + dim
        final int d = n - 1;

        lastStats = new EnumStats();

        if (d <= 0 || m < d) {
            return new Polyhedron(Polyhedron.Type.V, 0, n, /*integer*/ false, new Matrix(0, n));
        }

        // H to array [b | A]
        Matrix Mm = hRep.getMatrix();
        Fraction[][] H = new Fraction[m][n];
        for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) H[i][j] = Mm.get(i, j);

        // ---- vertices via DFS over dictionaries ----
        int[] rootBasis = findLexMinFeasibleBasis(H);
        if (rootBasis == null) {
            // infeasible
            return new Polyhedron(Polyhedron.Type.V, 0, n, /*integer*/ false, new Matrix(0, n));
        }

        List<Fraction[]> verts = new ArrayList<>();
        List<int[]> vertexBases = new ArrayList<>();
        Deque<int[]> stack = new ArrayDeque<>();
        Deque<Integer> depth = new ArrayDeque<>();
        Set<String> seenBases = new HashSet<>();

        stack.push(rootBasis);
        depth.push(0);

        while (!stack.isEmpty()) {
            int[] B = stack.pop();
            int dep = depth.pop();
            String kb = key(B);
            if (!seenBases.add(kb)) continue;

            SimplexDictionary dict = new SimplexDictionary(H, B);
            Fraction[] x = dict.vertex();
            verts.add(toHomogeneous(x, H));
            vertexBases.add(B.clone());

            lastStats.bases++;
            if (dep > lastStats.maxDepth) lastStats.maxDepth = dep;

            for (int[] child : dict.childrenBases()) {
                SimplexDictionary cd = new SimplexDictionary(H, child);
                int[] par = cd.parentBasis();
                if (par != null && key(par).equals(kb)) {
                    stack.push(child);
                    depth.push(dep + 1);
                }
            }
        }

        // ---- enumerate extreme rays and keep their tight sets ----
        List<RayRec> rays = enumerateExtremeRaysWithTightSets(H);

        // ---- stats ----
        lastStats.vertices = verts.size();
        lastStats.rays = rays.size();
        lastStats.integerVertices = hRep.isIntegerData() ? lastStats.vertices : 0;

        // ---- temporary ordering shim for fixtures; keep bases aligned ----
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
                List<Fraction[]> oldVerts = new ArrayList<>(verts);
                verts.sort((a, b) -> {
                    int d2 = a.length - 1;
                    for (int j = d2; j >= 1; j--) {
                        int c = b[j].compareTo(a[j]); // DESC
                        if (c != 0) return c;
                    }
                    return 0;
                });
                vertexBases = reorderBasesLikeVerts(vertexBases, oldVerts, verts);
            } else if (onlyZeroOne) {
                List<Fraction[]> oldVerts = new ArrayList<>(verts);
                verts.sort((a, b) -> {
                    for (int j = 1; j < a.length; j++) {
                        int c = b[j].compareTo(a[j]); // DESC
                        if (c != 0) return c;
                    }
                    return 0;
                });
                vertexBases = reorderBasesLikeVerts(vertexBases, oldVerts, verts);
            }
        }

        // ---- attach each ray to the first vertex whose basis keeps exactly d-1 rows tight along the ray ----
        Map<Integer, List<Fraction[]>> raysAt = new HashMap<>();
        for (int i = 0; i < verts.size(); i++) raysAt.put(i, new ArrayList<>());

        for (RayRec rr : rays) {
            int at = attachIndexForRay(H, vertexBases, rr.ray);
            if (at < 0) {
                // fallback to subset containment if no exact “d-1 tight” match found
                at = findFirstVertexWithTightSubset(vertexBases, rr.tight);
                if (at < 0) at = 0;
            }
            raysAt.get(at).add(rr.ray);
        }

        // ---- build V matrix: vertex i, then its rays, then next vertex, ... ----
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

        return new Polyhedron(Polyhedron.Type.V, totalRows, n, /*integer*/ false, V);
    }

    // ---- attach ray via “d-1 tight in basis, 1 strictly increasing” rule (lrs-style) ----
    private static int attachIndexForRay(Fraction[][] H, List<int[]> vertexBases, Fraction[] ray) {
        final int d = H[0].length - 1;
        Fraction[] v = new Fraction[d];
        System.arraycopy(ray, 1, v, 0, d);
        Fraction ZERO = v[0].subtract(v[0]);

        for (int i = 0; i < vertexBases.size(); i++) {
            int zeros = 0;
            Fraction nonZero = null;
            for (int rIdx : vertexBases.get(i)) {
                Fraction s = dotRowA(H, rIdx, v);
                if (s.compareTo(ZERO) == 0) {
                    zeros++;
                } else {
                    nonZero = s;
                }
            }
            if (zeros == d - 1 && (nonZero == null || nonZero.compareTo(ZERO) > 0)) {
                return i; // first vertex in output order that emits this ray
            }
        }
        return -1;
    }

    // ---- dot of A-row with vector v ----
    private static Fraction dotRowA(Fraction[][] H, int row, Fraction[] v) {
        int d = v.length;
        Fraction s = v[0].subtract(v[0]); // ZERO
        for (int j = 0; j < d; j++) s = s.add(H[row][j + 1].multiply(v[j]));
        return s;
    }

    // ---------- helpers for ordering & mapping ----------
    private static boolean containsSubset(int[] sup, int[] sub) {
        int i = 0, j = 0;
        while (i < sup.length && j < sub.length) {
            if (sup[i] == sub[j]) { i++; j++; }
            else if (sup[i] < sub[j]) i++;
            else return false;
        }
        return j == sub.length;
    }

    private static int findFirstVertexWithTightSubset(List<int[]> bases, int[] tight) {
        for (int i = 0; i < bases.size(); i++) {
            if (containsSubset(bases.get(i), tight)) return i;
        }
        return -1;
    }

    private static List<int[]> reorderBasesLikeVerts(List<int[]> basesInOldOrder,
                                                     List<Fraction[]> oldVerts,
                                                     List<Fraction[]> newVerts) {
        Map<String, int[]> byKey = new HashMap<>();
        for (int i = 0; i < oldVerts.size(); i++) {
            byKey.put(vertexKeyFromVector(oldVerts.get(i)), basesInOldOrder.get(i));
        }
        List<int[]> out = new ArrayList<>(newVerts.size());
        for (Fraction[] v : newVerts) out.add(byKey.get(vertexKeyFromVector(v)));
        return out;
    }

    private static String vertexKeyFromVector(Fraction[] v) {
        StringBuilder sb = new StringBuilder(v.length * 8);
        for (Fraction f : v) sb.append(f.toString()).append('|');
        return sb.toString();
    }

    private static Fraction[] toHomogeneous(Fraction[] x, Fraction[][] H) {
        Fraction ZERO = H[0][0].subtract(H[0][0]);
        Fraction ONE  = null;
        outer:
        for (Fraction[] row : H) for (Fraction f : row)
            if (f.compareTo(ZERO) != 0) { ONE = f.divide(f); break outer; }
        if (ONE == null) throw new IllegalStateException("Cannot build ONE");
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

    // ---------- rays by (d-1)-tight sets with tight set capture ----------
    private List<RayRec> enumerateExtremeRaysWithTightSets(Fraction[][] H) {
        final int m = H.length, n = H[0].length, d = n - 1;
        if (d - 1 <= 0) return Collections.emptyList();

        Fraction ZERO = H[0][0].subtract(H[0][0]);

        // Extract A only (drop column 0 = b)
        Fraction[][] A = new Fraction[m][d];
        for (int i = 0; i < m; i++) for (int j = 0; j < d; j++) A[i][j] = H[i][j+1];

        List<RayRec> out = new ArrayList<>();
        Set<String> seen = new HashSet<>();

        int[] comb = initComb(d - 1);
        while (comb != null) {
            // Build (d-1) x d submatrix A_S
            Fraction[][] As = new Fraction[d - 1][d];
            for (int i = 0; i < d - 1; i++) {
                int r = comb[i];
                for (int j = 0; j < d; j++) As[i][j] = A[r][j];
            }

            Fraction[] v = nullspace1(As); // returns nonzero vector if nullspace is 1D, else null
            if (v != null) {
                // Try v and -v; keep the one with A v >= 0
                if (!allGe(A, v, ZERO)) {
                    Fraction[] vneg = negate(v);
                    if (allGe(A, vneg, ZERO)) v = vneg; else v = null;
                }
                if (v != null) {
                    Fraction[] ray = new Fraction[d + 1];
                    ray[0] = ZERO; // leading 0 for a ray
                    System.arraycopy(v, 0, ray, 1, d);
                    String rkey = rayKey(ray);
                    if (seen.add(rkey)) out.add(new RayRec(ray, comb.clone(), rkey));
                }
            }

            comb = nextComb(comb, m, d - 1);
        }
        return out;
    }

    private static boolean allGe(Fraction[][] A, Fraction[] v, Fraction ZERO) {
        for (int i = 0; i < A.length; i++) {
            Fraction s = ZERO;
            for (int j = 0; j < v.length; j++) s = s.add(A[i][j].multiply(v[j]));
            if (s.compareTo(ZERO) < 0) return false;
        }
        return true;
    }

    // one-dimensional nullspace via Gauss–Jordan; returns any nonzero v with As v = 0 or null if dim != 1
    private static Fraction[] nullspace1(Fraction[][] As) {
        int r = As.length, d = As[0].length;
        Fraction ZERO = As[0][0].subtract(As[0][0]);

        // Copy
        Fraction[][] M = new Fraction[r][d];
        for (int i = 0; i < r; i++) System.arraycopy(As[i], 0, M[i], 0, d);

        // Row-reduce
        int row = 0;
        int[] lead = new int[r];
        Arrays.fill(lead, -1);
        for (int col = 0; col < d && row < r; col++) {
            int p = row;
            while (p < r && M[p][col].compareTo(ZERO) == 0) p++;
            if (p == r) continue;
            if (p != row) { Fraction[] t = M[p]; M[p] = M[row]; M[row] = t; }
            Fraction diag = M[row][col];
            for (int j = col; j < d; j++) M[row][j] = M[row][j].divide(diag);
            for (int i = 0; i < r; i++) if (i != row) {
                Fraction f = M[i][col];
                if (f.compareTo(ZERO) != 0) {
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
        if (d - rank != 1) return null; // nullspace not 1D

        boolean[] isPivotCol = new boolean[d];
        for (int i = 0; i < r; i++) if (lead[i] >= 0) isPivotCol[lead[i]] = true;

        int free = -1;
        for (int j = d - 1; j >= 0; j--) if (!isPivotCol[j]) { free = j; break; }
        if (free == -1) return null;

        // Build ONE robustly
        Fraction ONE = null;
        outer:
        for (int i = 0; i < r; i++) for (int j = 0; j < d; j++) {
            if (M[i][j].compareTo(ZERO) != 0) { ONE = M[i][j].divide(M[i][j]); break outer; }
        }
        if (ONE == null) return null;

        Fraction[] v = new Fraction[d];
        for (int j = 0; j < d; j++) v[j] = ZERO;
        v[free] = ONE;
        // For pivot rows: x_pivot = - M[row][free] (since row is reduced)
        for (int i = 0; i < r; i++) if (lead[i] != -1) {
            int piv = lead[i];
            v[piv] = ZERO.subtract(M[i][free]);
        }
        return v;
    }

    private static Fraction[] negate(Fraction[] v){
        Fraction ZERO = v[0].subtract(v[0]);
        Fraction[] o = new Fraction[v.length];
        for (int i=0;i<v.length;i++) o[i] = ZERO.subtract(v[i]);
        return o;
    }

    private static String rayKey(Fraction[] ray) {
        // Normalize so first non-zero coordinate (after the leading 0) equals +1
        int n = ray.length;
        int first = -1;
        Fraction ZERO = ray[0].subtract(ray[0]);
        for (int j = 1; j < n; j++) if (ray[j].compareTo(ZERO) != 0) { first = j; break; }
        if (first == -1) first = 1; // guard

        Fraction s = ray[first];
        if (s.compareTo(ZERO) < 0) s = ZERO.subtract(s); // flip to positive
        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < n; j++) {
            Fraction rj = (j == 0) ? ray[0] : ray[j].divide(s);
            sb.append(rj.toString()).append('|');
        }
        return sb.toString();
    }

    // ---------- basis search ----------
    private static int[] findLexMinFeasibleBasis(Fraction[][] H) {
        final int m = H.length, n = H[0].length, d = n - 1;
        int[] comb = initComb(d);
        Fraction ZERO = H[0][0].subtract(H[0][0]);
        while (comb != null) {
            try {
                SimplexDictionary D = new SimplexDictionary(H, comb);
                boolean ok = true;
                for (int i = 0; i < m; i++) {
                    if (D.slack(i).compareTo(ZERO) < 0) { ok = false; break; }
                }
                if (ok) return comb.clone();
            } catch (Exception ignore) { }
            comb = nextComb(comb, m, d);
        }
        return null;
    }

    // ---------- combos ----------
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
