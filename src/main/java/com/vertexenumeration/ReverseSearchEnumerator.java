package com.vertexenumeration;

import java.util.*;

/**
 * Reverse search enumerator (Java translation of lrslib/reverse.c).
 * Enumerates all vertices of an H-polyhedron via dictionary pivoting.
 */
final class ReverseSearchEnumerator {

    static final class Result {
        final List<Fraction[]> vertices;
        final List<Fraction[]> rays;
        final EnumStats stats;

        Result(List<Fraction[]> v, List<Fraction[]> r, EnumStats s) {
            vertices = v;
            rays = r;
            stats = s;
        }
    }


    private final Fraction[][] H; // inequalities: rows [b | A]
    private final int m, n, d;

    ReverseSearchEnumerator(Fraction[][] H) {
        this.H = H;
        this.m = H.length;
        this.n = (m == 0 ? 0 : H[0].length);
        this.d = Math.max(0, n - 1);
        if (m == 0 || n == 0 || d == 0)
            throw new IllegalArgumentException("Empty/degenerate matrix");
    }

    /** Run reverse search enumeration. */
    /** Run reverse search enumeration. */
    Result run() {
        EnumStats stats = new EnumStats();
        List<Fraction[]> verts = new ArrayList<>();
        Map<String, Fraction[]> uniqRays = new HashMap<>();

        int[] rootBasis = findFirstBasisPhaseI(H);
        if (rootBasis == null)
            return new Result(Collections.emptyList(), Collections.emptyList(), stats);

        Deque<int[]> stack = new ArrayDeque<>();
        Deque<Integer> depth = new ArrayDeque<>();
        Set<String> seen = new HashSet<>();

        stack.push(rootBasis);
        depth.push(0);

        while (!stack.isEmpty()) {
            int[] B = stack.pop();
            int dep = depth.pop();

            String kb = key(B);
            if (!seen.add(kb)) continue;

            SimplexDictionary dict = new SimplexDictionary(H, B);

            // ---- Vertex ----
            Fraction[] x = dict.vertex();
            Fraction[] homog = toHomogeneous(x);
            verts.add(homog);

            stats.bases++;
            stats.vertices++;
            if (isIntegerVertex(homog)) stats.integerVertices++;
            if (dep > stats.maxDepth) stats.maxDepth = dep;

            // ---- Rays ----
            for (Fraction[] ray : dict.rayDirections()) {
                Fraction[] canon = SimplexDictionary.canonicalizeRay(ray); // normalize
                String rk = Arrays.toString(canon); // safe key after canonicalization
                if (!uniqRays.containsKey(rk)) {
                    uniqRays.put(rk, canon);
                    stats.rays++;
                }
            }

            // ---- Children ----
            List<int[]> children = dict.childrenBases();
            // Reverse order before pushing → ensures DFS visits in lex order
            Collections.reverse(children);

            for (int[] child : children) {
                SimplexDictionary cd = new SimplexDictionary(H, child);
                int[] par = cd.parentBasis();
                if (par != null && key(par).equals(kb)) {
                    stack.push(child);
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



    // ---------- Utilities ----------

    /** Homogenize x → [1, x...] */
    private static Fraction[] toHomogeneous(Fraction[] x) {
        Fraction[] v = new Fraction[x.length + 1];
        v[0] = Fraction.ONE;
        System.arraycopy(x, 0, v, 1, x.length);
        return v;
    }

    /** Check if homogeneous vertex is integer. */
    private static boolean isIntegerVertex(Fraction[] v) {
        for (int j = 1; j < v.length; j++) {
            if (!v[j].denominator().equals(java.math.BigInteger.ONE))
                return false;
        }
        return true;
    }


    /** Phase I simplex: find a feasible basis like lrs_getfirstbasis. */
    /** Phase I simplex: find a feasible basis like lrs_getfirstbasis. */
    private static int[] findFirstBasisPhaseI(Fraction[][] H) {
        final int m = H.length, n = H[0].length, d = n - 1;

        // --- Step 1: Try trivial lex-min basis (first d rows) ---
        int[] basis = initComb(d);
        try {
            SimplexDictionary dict = new SimplexDictionary(H, basis);
            boolean feasible = true;
            for (int i = 0; i < m; i++) {
                if (dict.slack(i).compareTo(Fraction.ZERO) < 0) {
                    feasible = false;
                    break;
                }
            }
            if (feasible) {
                Arrays.sort(basis);
                return basis;
            }
        } catch (Exception ignore) {}

        // --- Step 2: Build Phase I system (add d artificial rows) ---
        Fraction[][] Hphase = new Fraction[m + d][n];
        for (int i = 0; i < m; i++) {
            Hphase[i] = Arrays.copyOf(H[i], n);
        }
        for (int i = 0; i < d; i++) {
            Hphase[m + i] = new Fraction[n];
            Arrays.fill(Hphase[m + i], Fraction.ZERO);
            Hphase[m + i][0] = Fraction.ONE;       // slack constant
            Hphase[m + i][i + 1] = Fraction.ONE;   // artificial identity
        }
        int mPhase = Hphase.length;

        // --- Step 3: Start with artificial basis (last d rows) ---
        int[] artBasis = new int[d];
        for (int i = 0; i < d; i++) artBasis[i] = m + i;

        SimplexDictionary dict;
        try {
            dict = new SimplexDictionary(Hphase, artBasis);
        } catch (RuntimeException ex) {
            System.err.println("*unrecoverable error: no nonsingular artificial basis");
            return null;
        }

        // --- Step 4: Pivot artificials out deterministically ---
        boolean progress = true;
        while (progress) {
            progress = false;

            // Check feasibility wrt original constraints
            boolean feasible = true;
            for (int i = 0; i < m; i++) {
                if (dict.slack(i).compareTo(Fraction.ZERO) < 0) {
                    feasible = false;
                    break;
                }
            }
            if (feasible) break;

            // Try to replace artificials with original constraints
            for (int e = 0; e < m; e++) {
                if (dict.slack(e).compareTo(Fraction.ZERO) < 0) {
                    int leave = dict.leavingFor(e);  // lex ratio rule
                    if (leave >= 0) {
                        int[] newBasis = dict.basis().clone();
                        newBasis[leave] = e;
                        Arrays.sort(newBasis);
                        dict = new SimplexDictionary(Hphase, newBasis);
                        progress = true;
                        break;
                    }
                }
            }

            if (!progress) {
                System.err.println("*unrecoverable error: infeasible system");
                return null;
            }
        }

        // --- Step 5: Drop artificial rows from basis ---
        int[] finalBasis = new int[d];
        int idx = 0;
        for (int b : dict.basis()) {
            if (b < m) { // keep only original rows
                if (idx < d) finalBasis[idx++] = b;
            }
        }
        if (idx != d) {
            System.err.println("*unrecoverable error: artificial still in basis");
            return null;
        }

        Arrays.sort(finalBasis);
        return finalBasis;
    }


    /**
     * Pick entering variable for Phase I: artificial column n.
     * Returns -1 if none can improve objective.
     */
    private static int chooseEnteringArtificial(SimplexDictionary dict, int artificialCol) {
        // Try to force artificial variable out of basis
        for (int e = 0; e < dict.basis().length; e++) {
            if (dict.basis()[e] == artificialCol) {
                return artificialCol;
            }
        }
        return -1; // artificial already out of basis
    }


    /** Ratio rule to pick leaving row. */
    private static int chooseLeaving(SimplexDictionary dict, int entering) {
        int leave = -1;
        Fraction bestT = null;
        int bestIdx = Integer.MAX_VALUE;
        Fraction ZERO = Fraction.ZERO;

        Fraction se = dict.slack(entering);
        if (se.compareTo(ZERO) >= 0) return -1;

        for (int i = 0; i < dict.basis().length; i++) {
            int r = dict.basis()[i];
            Fraction sr = dict.slack(r);
            if (sr.compareTo(ZERO) > 0) {
                Fraction ratio = sr.divide(se.abs());
                if (bestT == null || ratio.compareTo(bestT) < 0 ||
                        (ratio.compareTo(bestT) == 0 && i < bestIdx)) {
                    leave = i;
                    bestT = ratio;
                    bestIdx = i;
                }
            }
        }
        return leave;
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

    private static int[] nextComb(int[] a, int n, int k) {
        int i = k - 1;
        while (i >= 0 && a[i] == n - k + i) i--;
        if (i < 0) return null;
        a[i]++;
        for (int j = i + 1; j < k; j++) a[j] = a[j - 1] + 1;
        return a;
    }

    private static int lexCompareHomog(Fraction[] a, Fraction[] b) {
        for (int j = a.length - 1; j >= 1; j--) { // compare from last coordinate
            int cmp = a[j].compareTo(b[j]);
            if (cmp != 0) return -cmp; // lrslib uses descending
        }
        return 0;
    }


}
