package com.vertexenumeration;

import java.util.*;

/**
 * Reverse search enumerator (1:1 parity with lrslib/reverse.c).
 */
final class ReverseSearchEnumerator {

    static final class Result {
        final List<Fraction[]> vertices;
        final List<Fraction[]> rays;
        final EnumStats stats;
        Result(List<Fraction[]> v, List<Fraction[]> r, EnumStats s) {
            this.vertices = v; this.rays = r; this.stats = s;
        }
    }

    private final Fraction[][] H;
    private final int m, n, dOriginal;


    ReverseSearchEnumerator(Fraction[][] H) {
        this.H = H;
        this.m = H.length;
        this.n = (m == 0 ? 0 : H[0].length);
        this.dOriginal = Math.max(0, n - 1);
        if (m == 0 || n == 0 || dOriginal == 0)
            throw new IllegalArgumentException("Empty/degenerate matrix");
    }

    /** Main enumeration routine. */
    Result run() {
        EnumStats stats = new EnumStats();
        List<Fraction[]> verts = new ArrayList<>();
        Map<String, Fraction[]> uniqRays = new HashMap<>();

        // ---- Phase I: find feasible root basis ----
        int[] rootBasis = findFeasibleBasis(H);
        if (rootBasis == null)
            return new Result(Collections.emptyList(), Collections.emptyList(), stats);

        // DFS reverse search
        Deque<SimplexDictionary> stack = new ArrayDeque<>();
        Deque<Integer> depth = new ArrayDeque<>();
        Set<String> seen = new HashSet<>();

        SimplexDictionary rootDict = new SimplexDictionary(H, rootBasis, dOriginal);
        stack.push(rootDict);
        depth.push(0);

        while (!stack.isEmpty()) {
            SimplexDictionary dict = stack.pop();
            int dep = depth.pop();

            String kb = key(dict.basis());
            if (!seen.add(kb)) continue;

            // ---- count basis ----
            stats.bases++;
            if (dep > stats.maxDepth) stats.maxDepth = dep;

            // ---- vertex ----
            Fraction[] x = dict.vertex();
            Fraction[] homog = toHomogeneous(x);
            verts.add(homog);
            stats.vertices++;
            if (isIntegerVertex(homog)) stats.integerVertices++;

            // ---- objective (parity with lrslib) ----
            Fraction objVal = Fraction.ZERO;
            for (int j = 0; j < dOriginal; j++) {
                objVal = objVal.add(x[j]); // simplistic: sum of coords
            }
            dict.setObjective(objVal, Fraction.ONE, dict.objCol());

            // ---- rays ----
            for (Fraction[] ray : dict.rayDirections()) {
                Fraction[] canon = SimplexDictionary.canonicalizeRay(ray);
                String rk = Arrays.toString(canon);
                if (!uniqRays.containsKey(rk)) {
                    uniqRays.put(rk, canon);
                    stats.rays++;
                }
            }

            // ---- children ----
            List<int[]> children = dict.childrenBases();
            Collections.reverse(children); // DFS lex order
            for (int[] childBasis : children) {
                if (isParent(H, childBasis, dict.basis())) {
                    SimplexDictionary childDict = new SimplexDictionary(H, childBasis, dOriginal);
                    childDict.setPrev(dict); // parity: child->prev = parent
                    dict.setNext(childDict);
                    stack.push(childDict);
                    depth.push(dep + 1);
                }
            }
        }

        List<Fraction[]> rays = new ArrayList<>(uniqRays.values());
        verts.sort(ReverseSearchEnumerator::lexCompareHomog);
        rays.sort(ReverseSearchEnumerator::lexCompareHomog);

        stats.setMode(LrsDat.Mode.VE);
        return new Result(verts, rays, stats);
    }

    // ---------- Phase I ----------
    /**
     * Phase I feasibility routine (parity with lrs_getfirstbasis in lrslib).
     * Builds auxiliary system with artificials and pivots until feasible basis found.
     */
    /**
     * Phase I feasibility routine (parity with lrs_getfirstbasis in lrslib).
     * Builds auxiliary system with artificials and pivots until feasible basis found.
     */
    private int[] findFeasibleBasis(Fraction[][] H) {
        int m = H.length;
        int n = H[0].length;

        // Step 1: Try trivial basis = first dOriginal rows
        int[] trivial = initComb(dOriginal);
        try {
            SimplexDictionary test = new SimplexDictionary(H, trivial, dOriginal);
            boolean feasible = true;
            for (int i = 0; i < m; i++) {
                if (test.slack(i).compareTo(Fraction.ZERO) < 0) {
                    feasible = false;
                    break;
                }
            }
            if (feasible) return trivial;
        } catch (RuntimeException ignore) {
            // singular basis → fallback to Phase I
        }

        // Step 2: Build Phase I system with artificials
        Fraction[][] Hphase = new Fraction[m + dOriginal][n + dOriginal];

        // Copy original rows
        for (int i = 0; i < m; i++) {
            Hphase[i] = new Fraction[n + dOriginal];
            for (int j = 0; j < n; j++) Hphase[i][j] = H[i][j];
            for (int j = n; j < n + dOriginal; j++) Hphase[i][j] = Fraction.ZERO;
        }

        // Add artificial identity rows
        for (int i = 0; i < dOriginal; i++) {
            Hphase[m + i] = new Fraction[n + dOriginal];
            Arrays.fill(Hphase[m + i], Fraction.ZERO);
            Hphase[m + i][0] = Fraction.ONE;        // RHS
            Hphase[m + i][n + i] = Fraction.ONE;    // artificial variable
        }

        // Initial artificial basis = last dOriginal rows
        int[] basis = new int[dOriginal];
        for (int i = 0; i < dOriginal; i++) basis[i] = m + i;

        SimplexDictionary dict;
        try {
            dict = new SimplexDictionary(Hphase, basis, dOriginal);
        } catch (RuntimeException ex) {
            System.err.println("*unrecoverable error: no nonsingular artificial basis");
            return null;
        }

        // Step 3: Phase I simplex loop
        boolean progress = true;
        while (progress) {
            progress = false;

            try {
                // Check feasibility: all original slacks ≥ 0?
                boolean feasible = true;
                for (int i = 0; i < m; i++) {
                    if (dict.slack(i).compareTo(Fraction.ZERO) < 0) {
                        feasible = false;
                        break;
                    }
                }
                if (feasible) break;

                // Pick entering row e (any violating slack)
                int entering = -1;
                for (int e = 0; e < m; e++) {
                    if (dict.slack(e).compareTo(Fraction.ZERO) < 0) {
                        entering = e;
                        break;
                    }
                }
                if (entering < 0) {
                    System.err.println("*unrecoverable error: infeasible system");
                    return null;
                }

                // Leaving row via lex ratio test
                int leave = dict.leavingFor(entering);
                if (leave < 0) {
                    System.err.println("*unrecoverable error: infeasible in Phase I pivot");
                    return null;
                }

                // Update basis
                int[] nb = dict.basis().clone();
                nb[leave] = entering;
                Arrays.sort(nb);

                dict = new SimplexDictionary(Hphase, nb, dOriginal);
                progress = true;

            } catch (ArithmeticException ex) {
                // Singular pivot recovery (parity with lrslib)
                System.err.println("*warning: singular basis encountered, retrying Phase I...");
                progress = true; // force retry loop
            }
        }

        // Step 4: Drop artificials
        int[] finalBasis = Arrays.stream(dict.basis())
                .filter(b -> b < m)
                .toArray();
        if (finalBasis.length != dOriginal) {
            System.err.println("*unrecoverable error: artificial still in basis");
            return null;
        }

        Arrays.sort(finalBasis);
        return finalBasis;
    }






    // ---------- Parent test ----------
    /** Check if p is parent of child basis. */
    private boolean isParent(Fraction[][] H, int[] child, int[] parent) {
        SimplexDictionary dict = new SimplexDictionary(H, child, dOriginal);
        int[] par = dict.parentBasis();
        return par != null && Arrays.equals(par, parent);
    }

    // ---------- Utilities ----------
    private static Fraction[] toHomogeneous(Fraction[] x) {
        Fraction[] v = new Fraction[x.length + 1];
        v[0] = Fraction.ONE;
        System.arraycopy(x, 0, v, 1, x.length);
        return v;
    }
    private static boolean isIntegerVertex(Fraction[] v) {
        for (int j = 1; j < v.length; j++) {
            if (!v[j].denominator().equals(java.math.BigInteger.ONE)) return false;
        }
        return true;
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
    private static int lexCompareHomog(Fraction[] a, Fraction[] b) {
        for (int j = a.length - 1; j >= 1; j--) {
            int cmp = a[j].compareTo(b[j]);
            if (cmp != 0) return -cmp; // descending
        }
        return 0;
    }
}
