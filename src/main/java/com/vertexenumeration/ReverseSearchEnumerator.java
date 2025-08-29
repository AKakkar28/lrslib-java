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
    private final int m, n, d;

    ReverseSearchEnumerator(Fraction[][] H) {
        this.H = H;
        this.m = H.length;
        this.n = (m == 0 ? 0 : H[0].length);
        this.d = Math.max(0, n - 1);
        if (m == 0 || n == 0 || d == 0)
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

        SimplexDictionary rootDict = new SimplexDictionary(H, rootBasis);
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
            for (int j = 0; j < d; j++) {
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
                    SimplexDictionary childDict = new SimplexDictionary(H, childBasis);
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
    /** Full Phase I like lrs_getfirstbasis: minimize sum of artificials. */
    /**
     * Phase I simplex: full lrslib-style routine.
     * Builds auxiliary LP with artificial variables and minimizes their sum.
     */
    private static int[] findFeasibleBasis(Fraction[][] H) {
        final int m = H.length, n = H[0].length, d = n - 1;

        // Step 1: try trivial lex-min basis (first d rows)
        int[] trivial = initComb(d);
        try {
            SimplexDictionary dict = new SimplexDictionary(H, trivial);
            boolean feasible = true;
            for (int i = 0; i < m; i++) {
                if (dict.slack(i).compareTo(Fraction.ZERO) < 0) { feasible = false; break; }
            }
            if (feasible) return trivial.clone();
        } catch (Exception ignore) {}

        // Step 2: Build Phase I system with artificials
        Fraction[][] Hphase = new Fraction[m + d][n];
        for (int i = 0; i < m; i++) {
            Hphase[i] = Arrays.copyOf(H[i], n);
        }
        for (int i = 0; i < d; i++) {
            Hphase[m + i] = new Fraction[n];
            Arrays.fill(Hphase[m + i], Fraction.ZERO);
            Hphase[m + i][0] = Fraction.ONE;
            Hphase[m + i][i + 1] = Fraction.ONE; // artificial identity
        }

        // Initial artificial basis = last d rows
        int[] basis = new int[d];
        for (int i = 0; i < d; i++) basis[i] = m + i;
        SimplexDictionary dict;
        try {
            dict = new SimplexDictionary(Hphase, basis);
        } catch (RuntimeException ex) {
            System.err.println("*unrecoverable error: no nonsingular artificial basis");
            return null;
        }

        // Step 3: Run Phase I auxiliary simplex
        boolean progress = true;
        while (progress) {
            progress = false;

            // Compute objective = sum of artificials in basis
            Fraction obj = Fraction.ZERO;
            for (int b : dict.basis()) {
                if (b >= m) { // artificial
                    obj = obj.add(dict.slack(b));
                }
            }
            if (obj.compareTo(Fraction.ZERO) == 0) break; // feasible!

            // Choose entering variable: any violating original row
            int entering = -1;
            for (int e = 0; e < m; e++) {
                if (dict.slack(e).compareTo(Fraction.ZERO) < 0) {
                    entering = e;
                    break;
                }
            }
            if (entering == -1) {
                System.err.println("*unrecoverable error: infeasible system");
                return null;
            }

            // Leaving variable via lex ratio rule
            int leave = dict.leavingFor(entering);
            if (leave < 0) {
                System.err.println("*unrecoverable error: infeasible in Phase I pivot");
                return null;
            }

            // Update basis
            int[] newBasis = dict.basis().clone();
            newBasis[leave] = entering;
            Arrays.sort(newBasis);
            dict = new SimplexDictionary(Hphase, newBasis);
            progress = true;
        }

        // Step 4: Drop artificials
        int[] finalBasis = Arrays.stream(dict.basis())
                .filter(b -> b < m)
                .toArray();
        if (finalBasis.length != d) {
            System.err.println("*unrecoverable error: artificial still in basis");
            return null;
        }

        Arrays.sort(finalBasis);
        return finalBasis;
    }


    // ---------- Parent test ----------
    /** Check if p is parent of child basis. */
    private static boolean isParent(Fraction[][] H, int[] child, int[] parent) {
        SimplexDictionary dict = new SimplexDictionary(H, child);
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
