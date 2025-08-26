package com.vertexenumeration;

import java.util.*;

/** Exact simplex dictionary for H: b + A x >= 0 (no Phase-I here yet).
 *  Basis = d tight rows. Computes vertex, slacks, and lexicographic pivots.
 */
final class SimplexDictionary {
    private final Fraction[][] H; // m x n (n=1+d), columns [b | A]
    private final int m, n, d;

    private int[] basis;          // size d, sorted row indices
    private Fraction[][] Binv;    // d x d exact inverse of A_B
    private Fraction[] x;         // current vertex, length d

    SimplexDictionary(Fraction[][] H, int[] basisRows) {
        this.H = H;
        this.m = H.length;
        this.n = H[0].length;
        this.d = n - 1;
        this.basis = basisRows.clone();
        Arrays.sort(this.basis);
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
    }

    /** Slack of row i at current x: b_i + a_i·x */
    Fraction slack(int i) {
        Fraction s = H[i][0];
        for (int j = 0; j < d; j++) s = s.add(H[i][j+1].multiply(x[j]));
        return s;
    }

    /** Children in lexicographic order by lrs's rule:
     *  For each entering e, choose leaving via lex ratio rule. */
    // Children in lexicographic order using correct edge construction.
// For each nonbasic row e and each basic position l:
//   u := column l of B^{-1}   (i.e., B u = e_l)
//   denom := a_e · u
//   require denom < 0  (so a_l·dx = -1/denom > 0)
//   step t* = s_e  (because a_e·dx = -1)
//   check feasibility for all nonbasics j ≠ e: s_j + t* * (a_j·dx) >= 0
// If feasible, neighbor basis = (basis with row at position l replaced by e)
    List<int[]> childrenBases() {
        List<int[]> out = new ArrayList<>();
        BitSet inB = inBasis();
        Fraction ZERO = H[0][0].subtract(H[0][0]);

        for (int e = 0; e < m; e++) if (!inB.get(e)) {
            Fraction se = slack(e);               // current slack of entering row
            for (int l = 0; l < d; l++) {
                Fraction[] u = columnOfBinv(l);   // B u = e_l
                Fraction denom = dotRowA(e, u);   // a_e · u
                if (denom.compareTo(ZERO) >= 0) continue;     // need denom < 0
                boolean ok = true;

                // a_j·dx = -(a_j·u)/denom ; step t* = se
                for (int j = 0; j < m && ok; j++) {
                    if (inB.get(j)) continue;     // nonbasics only
                    if (j == e) continue;         // entering becomes tight
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


    /** Lex parent: lexicographically smallest neighbor strictly < this basis. */
    // Parent is the lexicographically smallest neighbor strictly less than this basis
    int[] parentBasis() {
        int[] best = null;
        for (int[] nb : childrenBases()) {
            if (lexCompare(nb, basis) < 0) {
                if (best == null || lexCompare(nb, best) < 0) best = nb;
            }
        }
        return best;
    }


    /** Compute leaving row index in basis for entering e using a lexicographic ratio rule. */
    private int leavingFor(int e) {
        // Direction: B * dx = -a_e
        Fraction[] a_e = new Fraction[d];
        for (int j = 0; j < d; j++) a_e[j] = H[e][j + 1];
        Fraction[] dx = solveWithBinv(negate(a_e));

        int leave = -1;
        Fraction bestT = null;
        int bestIdx = Integer.MAX_VALUE;
        Fraction ZERO = zero(a_e[0]);

        // For each basic row, compute a_r · dx (change in its slack).
        // Choose the lexicographically best ratio (t, idx) where t = se / (a_e · dx)
        // and a_r · dx > 0 (keeps feasibility when moving to make s_e -> 0).
        for (int i = 0; i < d; i++) {
            int r = basis[i];
            Fraction delta = dotRowA(r, dx);              // a_r · dx
            if (delta.compareTo(ZERO) > 0) {
                Fraction se = slack(e);                   // current nonbasic slack
                Fraction denom = dotRowA(e, dx);         // a_e · dx
                if (denom.compareTo(ZERO) == 0) continue; // parallel, no pivot
                Fraction t = se.divide(denom);           // step to drive s_e to 0
                if (betterRatio(t, i, bestT, bestIdx)) {
                    leave = i;
                    bestT = t;
                    bestIdx = i;
                }
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


    private static Fraction[] solve(Fraction[][] A, Fraction[] b){
        int n = b.length;
        Fraction ZERO = b[0].subtract(b[0]);
        Fraction[][] M = new Fraction[n][n+1];
        for (int i=0;i<n;i++){ System.arraycopy(A[i],0,M[i],0,n); M[i][n]=b[i]; }
        int r = 0;
        for (int c = 0; c < n && r < n; c++){
            int p=r; while(p<n && M[p][c].compareTo(ZERO)==0) p++;
            if (p==n) continue;
            if (p!=r){ Fraction[] t=M[p]; M[p]=M[r]; M[r]=t; }
            Fraction diag = M[r][c];
            for (int j=c;j<=n;j++) M[r][j]=M[r][j].divide(diag);
            for (int i2=0;i2<n;i2++) if (i2!=r){
                Fraction f = M[i2][c];
                if (f.compareTo(ZERO)!=0) for (int j=c;j<=n;j++) M[i2][j]=M[i2][j].subtract(f.multiply(M[r][j]));
            }
            r++;
        }
        Fraction[] x = new Fraction[n];
        for (int i=0;i<n;i++){ int lead=-1; for(int j=0;j<n;j++) if(M[i][j].compareTo(ZERO)!=0){lead=j;break;} if(lead==-1) return null; x[lead]=M[i][n]; }
        for (int i=0;i<n;i++) if (x[i]==null) return null;
        return x;
    }

    private static Fraction[][] invert(Fraction[][] B){
        int n = B.length;
        // Build ZERO and ONE robustly
        Fraction ZERO = B[0][0].subtract(B[0][0]);
        Fraction ONE = null;
        outer:
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (B[i][j].compareTo(ZERO) != 0) {
                    ONE = B[i][j].divide(B[i][j]); // non-zero / itself = 1
                    break outer;
                }
            }
        }
        if (ONE == null) {
            // Matrix is all zeros -> singular; let caller handle via later failure
            // but avoid 0/0 here by fabricating ONE as ZERO-ZERO+ZERO... can't.
            // Better: throw to signal singular immediately.
            throw new ArithmeticException("Singular matrix in invert(): all zeros");
        }

        // Augment [B | I]
        Fraction[][] M = new Fraction[n][2*n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) M[i][j] = B[i][j];
            for (int j = 0; j < n; j++) M[i][n + j] = (i == j) ? ONE : ZERO;
        }

        // Gauss–Jordan
        int r = 0;
        for (int c = 0; c < n && r < n; c++) {
            int p = r;
            while (p < n && M[p][c].compareTo(ZERO) == 0) p++;
            if (p == n) continue;             // no pivot in this column
            if (p != r) { Fraction[] t = M[p]; M[p] = M[r]; M[r] = t; }

            Fraction diag = M[r][c];
            // normalize pivot row
            for (int j = c; j < 2*n; j++) M[r][j] = M[r][j].divide(diag);

            // eliminate others
            for (int i2 = 0; i2 < n; i2++) if (i2 != r) {
                Fraction f = M[i2][c];
                if (f.compareTo(ZERO) != 0) {
                    for (int j = c; j < 2*n; j++) {
                        M[i2][j] = M[i2][j].subtract(f.multiply(M[r][j]));
                    }
                }
            }
            r++;
        }

        // Extract inverse
        Fraction[][] inv = new Fraction[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                inv[i][j] = M[i][n + j];
        return inv;
    }


    private BitSet inBasis(){ BitSet bs=new BitSet(m); for(int r:basis) bs.set(r); return bs; }
    private static int lexCompare(int[] a, int[] b) {
        for (int i = 0; i < a.length; i++) {
            int c = Integer.compare(a[i], b[i]);
            if (c != 0) return c;
        }
        return 0;
    }
    private static int lexFrac(Fraction[] A, Fraction[] B){
        for (int i=0;i<A.length;i++){ int c = A[i].compareTo(B[i]); if (c!=0) return c; }
        return 0;
    }
}
