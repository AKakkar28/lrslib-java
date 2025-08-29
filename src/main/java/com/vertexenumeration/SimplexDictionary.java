package com.vertexenumeration;

import java.util.*;

/** Exact simplex dictionary for H: b + A x >= 0 (no Phase-I here yet).
 *  Basis = d tight rows. Computes vertex, slacks, and lexicographic pivots.
 */
final class SimplexDictionary {
    // Tableau: H is m x n, with n = d+1 (RHS | A)
    private final Fraction[][] H;
    private final int m, n, d;

    // --- Parity fields from lrs_dic_struct ---
    private final int[] Row;       // row permutations
    private final int[] Col;       // col permutations
    private final int[] B;         // explicit basis bookkeeping
    private final int[] C;         // explicit cobasis bookkeeping
    private final boolean[] lexFlag; // parity: used for lexmin pivot choice

    private int depth;             // recursion depth
    private int pivotRow, pivotCol;// last pivot row/col
    private int lastDecisionVar;   // last entering variable

    // Objective function (lrslib stores numerator/denominator explicitly)
    private Fraction objNum;
    private Fraction objDen;
    private int objCol;            // objective column index

    private Dictionary.LPStatus status; // dictionary feasibility status

    // Reverse search tree links (lrslib uses prev/next pointers)
    private SimplexDictionary prev;
    private SimplexDictionary next;

    // Current basis & factorization
    private int[] basis;        // size d, sorted row indices
    private Fraction[][] Binv;  // d x d exact inverse of A_B
    private Fraction[] x;       // current vertex, length d


    SimplexDictionary(Fraction[][] H, int[] basisVars) {
        this.H = H;
        this.m = H.length;
        this.n = H[0].length;
        this.d = n - 1;

        // Basis/cobasis setup
        this.basis = basisVars.clone();
        Arrays.sort(this.basis);

        List<Integer> cob = new ArrayList<>();
        for (int j = 0; j < n; j++) {
            if (Arrays.binarySearch(basis, j) < 0) {
                cob.add(j);
            }
        }
        this.C = cob.stream().mapToInt(Integer::intValue).toArray();
        this.B = Arrays.copyOf(basis, basis.length);

        // Permutations
        this.Row = new int[m];
        this.Col = new int[n];
        for (int i = 0; i < m; i++) Row[i] = i;
        for (int j = 0; j < n; j++) Col[j] = j;

        // Lex flags
        this.lexFlag = new boolean[n];

        // Init values
        this.depth = 0;
        this.pivotRow = -1;
        this.pivotCol = -1;
        this.lastDecisionVar = -1;
        this.objCol = 0;
        this.objNum = Fraction.ZERO;
        this.objDen = Fraction.ONE;
        this.status = Dictionary.LPStatus.RUNNING;

        // Reverse search links null initially
        this.prev = null;
        this.next = null;

        refactor();
    }


    int[] basis() { return basis.clone(); }
    Fraction[] vertex() { return x.clone(); }

    /** Recompute B^{-1} and x from basis. */
    private void refactor() {
        // Build B and b_B
        Fraction[][] B = new Fraction[d][d];
        Fraction[] bB = new Fraction[d];
        for (int i = 0; i < d; i++) {
            int r = basis[i];
            bB[i] = H[r][0];
            for (int j = 0; j < d; j++) B[i][j] = H[r][j+1];
        }
        // Binv by solving B * Binv = I
        this.Binv = invert(B);
        // x solves B x = -b_B
        this.x = solve(B, negate(bB));
        if (this.x == null) {
            throw new RuntimeException("Refactor failed: singular basis " + Arrays.toString(basis));
        }

    }

    /** Slack of row i at current x: b_i + a_i·x */
    Fraction slack(int i) {
        Fraction s = H[i][0];
        for (int j = 0; j < d; j++) s = s.add(H[i][j+1].multiply(x[j]));
        return s;
    }

    public int depth() { return depth; }
    public int pivotRow() { return pivotRow; }
    public int pivotCol() { return pivotCol; }
    public int lastDecisionVar() { return lastDecisionVar; }
    public int objCol() { return objCol; }
    public Fraction objNum() { return objNum; }
    public Fraction objDen() { return objDen; }
    public Dictionary.LPStatus status() { return status; }
    public int[] rowPerm() { return Arrays.copyOf(Row, Row.length); }
    public int[] colPerm() { return Arrays.copyOf(Col, Col.length); }
    public boolean[] lexFlags() { return Arrays.copyOf(lexFlag, lexFlag.length); }
    public SimplexDictionary getPrev() { return prev; }
    public SimplexDictionary getNext() { return next; }
    public void setPrev(SimplexDictionary p) { this.prev = p; }
    public void setNext(SimplexDictionary n) { this.next = n; }

    // --- Methods for objective parity ---
    public void setObjective(Fraction num, Fraction den, int colIndex) {
        this.objNum = num;
        this.objDen = den;
        this.objCol = colIndex;
    }

    /** Children in lexicographic order by lrs's rule.
     * For each entering row e, choose leaving via lex ratio rule. */
    List<int[]> childrenBases() {
        List<int[]> out = new ArrayList<>();
        Fraction ZERO = H[0][0].subtract(H[0][0]);

        for (int e = 0; e < m; e++) {
            // skip if already in basis
            boolean inB = false;
            for (int b : basis) if (b == e) { inB = true; break; }
            if (inB) continue;

            Fraction se = slack(e);
            for (int l = 0; l < d; l++) {
                Fraction[] u = columnOfBinv(l);
                Fraction denom = dotRowA(e, u);
                if (denom.compareTo(ZERO) >= 0) continue;

                boolean ok = true;
                for (int j = 0; j < m && ok; j++) {
                    if (j == e) continue;
                    boolean inBasis = false;
                    for (int b : basis) if (b == j) { inBasis = true; break; }
                    if (inBasis) continue;

                    Fraction aj_u = dotRowA(j, u);
                    Fraction ajdx = ZERO.subtract(aj_u).divide(denom);
                    if (slack(j).add(se.multiply(ajdx)).compareTo(ZERO) < 0) ok = false;
                }

                if (ok) {
                    int[] nb = basis.clone();
                    nb[l] = e;
                    Arrays.sort(nb);
                    out.add(nb);
                }
            }
        }

        out.sort(SimplexDictionary::lexCompare);
        return out;
    }


    /** Lex parent = lexicographically smallest neighbor strictly less than this basis. */
    int[] parentBasis() {
        int[] best = null;
        for (int[] nb : childrenBases()) {
            if (lexCompare(nb, basis) < 0) {
                if (best == null || lexCompare(nb, best) < 0) {
                    best = nb;
                }
            }
        }
        return best;
    }


    /** Slack for a given nonbasic variable (column index). */
    private Fraction slackVar(int varCol) {
        // RHS + column contribution
        Fraction s = H[0][0].subtract(H[0][0]); // ZERO
        for (int i = 0; i < d; i++) {
            s = s.add(H[basis[i]][varCol].multiply(x[i]));
        }
        return s;
    }

    /** Dot product of a column (variable index) with a direction vector u. */
    private Fraction colDot(int varCol, Fraction[] u) {
        Fraction s = H[0][0].subtract(H[0][0]); // ZERO
        for (int i = 0; i < d; i++) {
            s = s.add(H[basis[i]][varCol].multiply(u[i]));
        }
        return s;
    }





    /**
     * Lexicographic ratio test (lrslib’s lrs_lexminratio).
     * Returns index of leaving row in the basis, or -1 if none.
     */
    public int leavingFor(int e) {
        Fraction[] a_e = new Fraction[d];
        for (int j = 0; j < d; j++) a_e[j] = H[e][j + 1];

        // Direction dx = -B^{-1} a_e
        Fraction[] dx = solveWithBinv(negate(a_e));
        Fraction se = slack(e);

        int leave = -1;
        Fraction[] bestLex = null;

        for (int i = 0; i < d; i++) {
            int r = basis[i];

            // Δs_r = a_r · dx
            Fraction delta = dotRowA(r, dx);
            if (delta.compareTo(Fraction.ZERO) <= 0) continue; // must improve slack

            // Lex ratio vector: [ t , row of B^{-1} ]
            // where t = s_e / delta
            Fraction t = se.divide(delta);
            Fraction[] ratioVec = new Fraction[d + 1];
            ratioVec[0] = t;
            System.arraycopy(Binv[i], 0, ratioVec, 1, d);

            if (bestLex == null || lexFrac(ratioVec, bestLex) < 0) {
                bestLex = ratioVec;
                leave = i;
            }
        }
        return leave;
    }
    private static boolean betterRatio(Fraction t, int idx, Fraction bestT, int bestIdx) {
        if (bestT == null) return true;
        int c = t.compareTo(bestT);
        if (c != 0) return c < 0;     // smaller step first
        return idx < bestIdx;         // deterministic tie-break without Fraction.valueOf
    }



    // ---------- math helpers (exact) ----------

    private static Fraction zero(Fraction f){ return f.subtract(f); }
    private static Fraction[] negate(Fraction[] v){ Fraction[] o=new Fraction[v.length]; for(int i=0;i<v.length;i++) o[i]=zero(v[i]).subtract(v[i]); return o; }

    private Fraction dotRowA(int r, Fraction[] v){
        Fraction s = zero(v[0]);
        for (int j = 0; j < d; j++) s = s.add(H[r][j+1].multiply(v[j]));
        return s;
    }

    private Fraction[] solveWithBinv(Fraction[] rhs){
        // dx = Binv * rhs
        Fraction[] out = new Fraction[d];
        for (int i = 0; i < d; i++) {
            Fraction s = zero(rhs[0]);
            for (int j = 0; j < d; j++) s = s.add(Binv[i][j].multiply(rhs[j]));
            out[i] = s;
        }
        return out;
    }

    private Fraction[] columnOfBinv(int col) {
        Fraction[] u = new Fraction[d];
        for (int i = 0; i < d; i++) u[i] = Binv[i][col];  // B^{-1} * e_col
        return u;
    }


    private static Fraction[] solve(Fraction[][] A, Fraction[] b) {
        int n = b.length;
        Fraction ZERO = b[0].subtract(b[0]);
        Fraction[][] M = new Fraction[n][n + 1];

        // build augmented matrix
        for (int i = 0; i < n; i++) {
            System.arraycopy(A[i], 0, M[i], 0, n);
            M[i][n] = b[i];
        }

        int[] colPerm = new int[n];
        for (int j = 0; j < n; j++) colPerm[j] = j;

        // Gaussian elimination with full pivoting
        for (int k = 0; k < n; k++) {
            // find pivot (max abs entry among submatrix)
            int pivRow = -1, pivCol = -1;
            outer:
            for (int i = k; i < n; i++) {
                for (int j = k; j < n; j++) {
                    if (M[i][j].compareTo(ZERO) != 0) {
                        pivRow = i; pivCol = j; break outer;
                    }
                }
            }
            if (pivRow == -1) throw new ArithmeticException("Singular system in solve()");

            // swap rows
            if (pivRow != k) {
                Fraction[] tmp = M[pivRow]; M[pivRow] = M[k]; M[k] = tmp;
            }
            // swap columns (track permutation)
            if (pivCol != k) {
                for (int i = 0; i < n; i++) {
                    Fraction tmp = M[i][pivCol]; M[i][pivCol] = M[i][k]; M[i][k] = tmp;
                }
                int tmpIdx = colPerm[pivCol]; colPerm[pivCol] = colPerm[k]; colPerm[k] = tmpIdx;
            }

            // normalize pivot row
            Fraction diag = M[k][k];
            for (int j = k; j <= n; j++) M[k][j] = M[k][j].divide(diag);

            // eliminate others
            for (int i = 0; i < n; i++) if (i != k) {
                Fraction f = M[i][k];
                if (f.compareTo(ZERO) != 0) {
                    for (int j = k; j <= n; j++) {
                        M[i][j] = M[i][j].subtract(f.multiply(M[k][j]));
                    }
                }
            }
        }

        // extract solution with column permutation reversed
        Fraction[] x = new Fraction[n];
        for (int i = 0; i < n; i++) x[colPerm[i]] = M[i][n];
        return x;
    }



    private static Fraction[][] invert(Fraction[][] B) {
        int n = B.length;
        Fraction ZERO = B[0][0].subtract(B[0][0]);

        // augmented [B | I]
        Fraction[][] M = new Fraction[n][2 * n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) M[i][j] = B[i][j];
            for (int j = 0; j < n; j++) M[i][n + j] = (i == j ? Fraction.ONE : ZERO);
        }

        int[] colPerm = new int[n];
        for (int j = 0; j < n; j++) colPerm[j] = j;

        // full pivoting
        for (int k = 0; k < n; k++) {
            int pivRow = -1, pivCol = -1;
            outer:
            for (int i = k; i < n; i++) {
                for (int j = k; j < n; j++) {
                    if (!M[i][j].equals(ZERO)) {
                        pivRow = i; pivCol = j; break outer;
                    }
                }
            }
            if (pivRow == -1) throw new ArithmeticException("Singular matrix in invert()");

            // row swap
            if (pivRow != k) {
                Fraction[] tmp = M[pivRow]; M[pivRow] = M[k]; M[k] = tmp;
            }
            // col swap (and track)
            if (pivCol != k) {
                for (int i = 0; i < n; i++) {
                    Fraction tmp = M[i][pivCol]; M[i][pivCol] = M[i][k]; M[i][k] = tmp;
                }
                int tmpIdx = colPerm[pivCol]; colPerm[pivCol] = colPerm[k]; colPerm[k] = tmpIdx;
            }

            // normalize pivot row
            Fraction diag = M[k][k];
            for (int j = k; j < 2 * n; j++) M[k][j] = M[k][j].divide(diag);

            // eliminate
            for (int i = 0; i < n; i++) if (i != k) {
                Fraction f = M[i][k];
                if (!f.equals(ZERO)) {
                    for (int j = k; j < 2 * n; j++) {
                        M[i][j] = M[i][j].subtract(f.multiply(M[k][j]));
                    }
                }
            }
        }

        // extract inverse with column permutation
        Fraction[][] inv = new Fraction[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                inv[colPerm[j]][i] = M[i][n + j]; // note permuted
            }
        }
        return inv;
    }



    private BitSet inBasis(){ BitSet bs=new BitSet(m); for(int r:basis) bs.set(r); return bs; }
    private static int lexCompare(int[] a, int[] b){ for(int i=0;i<a.length;i++){ int c=Integer.compare(a[i], b[i]); if(c!=0) return c; } return 0; }
    private static int lexFrac(Fraction[] A, Fraction[] B){
        for (int i=0;i<A.length;i++){ int c = A[i].compareTo(B[i]); if (c!=0) return c; }
        return 0;
    }

    /** Returns list of extreme rays from this basis (homogenized form: [0 | r]). */
    /** Returns list of extreme rays from this basis (homogenized form: [0 | r]). */
    List<Fraction[]> rayDirections() {
        List<Fraction[]> rays = new ArrayList<>();
        BitSet inB = inBasis();
        Fraction ZERO = H[0][0].subtract(H[0][0]);

        for (int e = 0; e < m; e++) if (!inB.get(e)) {
            // Candidate direction: solve B * dx = -a_e
            Fraction[] a_e = new Fraction[d];
            for (int j = 0; j < d; j++) a_e[j] = H[e][j + 1];
            Fraction[] dx = solveWithBinv(negate(a_e));

            // Homogenized ray: [0 | dx]
            Fraction[] ray = new Fraction[d + 1];
            ray[0] = ZERO;
            System.arraycopy(dx, 0, ray, 1, d);

            // Feasibility check:
            // - a_e·dx == 0 (entering row tight)
            // - all other nonbasic rows j ≠ e: a_j·dx >= 0
            boolean feasible = true;
            if (!dotRowA(e, dx).equals(ZERO)) feasible = false;
            for (int j = 0; j < m && feasible; j++) {
                if (inB.get(j) || j == e) continue;
                if (dotRowA(j, dx).compareTo(ZERO) < 0) feasible = false;
            }

            if (feasible) {
                rays.add(canonicalizeRay(ray));
            }
        }
        return rays;
    }


    /** Normalize ray: divide out gcd, flip sign if needed so first nonzero > 0. */
    /** Normalize ray: divide out gcd, flip sign so first nonzero > 0. */
    public static Fraction[] canonicalizeRay(Fraction[] r) {
        int n = r.length;
        Fraction ZERO = r[0].subtract(r[0]);

        // Find first nonzero coordinate (ignoring homogeneous 0)
        int first = -1;
        for (int j = 1; j < n; j++) {
            if (!r[j].equals(ZERO)) { first = j; break; }
        }
        if (first == -1) return r.clone(); // zero ray (shouldn't happen)

        // Collect numerators/denominators to compute gcd
        java.math.BigInteger g = java.math.BigInteger.ZERO;
        for (int j = 1; j < n; j++) {
            java.math.BigInteger num = r[j].numerator().abs();
            if (!num.equals(java.math.BigInteger.ZERO)) {
                g = (g.equals(java.math.BigInteger.ZERO)) ? num : g.gcd(num);
            }
        }

        Fraction[] out = new Fraction[n];
        if (!g.equals(java.math.BigInteger.ZERO)) {
            for (int j = 0; j < n; j++) {
                out[j] = new Fraction(r[j].numerator().divide(g), r[j].denominator());
            }
        } else {
            System.arraycopy(r, 0, out, 0, n);
        }

        // Flip sign if first nonzero < 0
        if (out[first].compareTo(ZERO) < 0) {
            for (int j = 0; j < n; j++) {
                out[j] = out[j].negate();
            }
        }

        return out;
    }



}
