package com.vertexenumeration;

import java.util.Arrays;
import java.util.Objects;

/**
 * Simplex dictionary (tableau) with lrs-style lexicographic ratio test.
 *
 * Indices:
 *   rows: 0 = objective, 1..m = constraints
 *   cols: 0 = RHS, 1..n = non-basic variable columns
 *
 * basis[r]   = var id basic in row r (r>=1)
 * cobasis[c] = var id non-basic in column c (c>=1)
 */
final class Dictionary implements LrsDic {

    private final int m, n;
    private final Fraction[][] T;   // (m+1) x (n+1)
    private final int[] basis;      // size m+1, basis[0] unused
    private final int[] cobasis;    // size n+1, cobasis[0] unused

    // Default to LEXICOGRAPHIC to match lrs; you can set BLAND for debugging.
    private PivotRule pivotRule = PivotRule.LEXICOGRAPHIC;

    private LPStatus status = LPStatus.RUNNING;

    Dictionary(int m, int n, Fraction[][] tableau, int[] basis, int[] cobasis) {
        this.m = m;
        this.n = n;
        this.T = deepCopy(tableau);
        this.basis = Arrays.copyOf(basis, basis.length);
        this.cobasis = Arrays.copyOf(cobasis, cobasis.length);
        sanity();
    }

    public static Dictionary of(Fraction[][] tableau, int[] basis, int[] cobasis) {
        int m = tableau.length - 1;
        int n = tableau[0].length - 1;
        return new Dictionary(m, n, tableau, basis, cobasis);
    }

    private void sanity() {
        if (T.length != m + 1) throw new IllegalArgumentException("tableau row count mismatch");
        for (Fraction[] row : T) {
            if (row.length != n + 1) throw new IllegalArgumentException("tableau col count mismatch");
            for (int j = 0; j <= n; j++) if (row[j] == null) throw new IllegalArgumentException("null Fraction");
        }
        if (basis.length != m + 1) throw new IllegalArgumentException("basis length mismatch");
        if (cobasis.length != n + 1) throw new IllegalArgumentException("cobasis length mismatch");
    }

    private static Fraction[][] deepCopy(Fraction[][] a) {
        Fraction[][] c = new Fraction[a.length][];
        for (int i = 0; i < a.length; i++) c[i] = Arrays.copyOf(a[i], a[i].length);
        return c;
    }

    // --- accessors
    public int m() { return m; }
    public int n() { return n; }
    public LPStatus status() { return status; }
    public int[] basis() { return Arrays.copyOf(basis, basis.length); }
    public int[] cobasis() { return Arrays.copyOf(cobasis, cobasis.length); }
    public Fraction[][] tableau() { return deepCopy(T); }
    public void setPivotRule(PivotRule rule) { this.pivotRule = Objects.requireNonNull(rule); }

    private static boolean isNeg(Fraction x)  { return x.compareTo(Fraction.ZERO) < 0; }
    private static boolean isPos(Fraction x)  { return x.compareTo(Fraction.ZERO) > 0; }
    private static boolean isZero(Fraction x) { return x.compareTo(Fraction.ZERO) == 0; }

    /** Feasible if RHS in all constraint rows is >= 0. */
    public boolean isFeasible() {
        for (int r = 1; r <= m; r++) if (isNeg(T[r][0])) return false;
        return true;
    }

    private boolean anyNegativeReducedCost() {
        for (int c = 1; c <= n; c++) if (isNeg(T[0][c])) return true;
        return false;
    }

    /** Single simplex iteration (minimization with reduced costs in row 0). */
    public void step() {
        if (!anyNegativeReducedCost()) { status = LPStatus.OPTIMAL; return; }

        int e = chooseEnteringColumn();
        if (e == -1) { status = LPStatus.OPTIMAL; return; }

        int r = chooseLeavingRow(e);
        if (r == -1) { status = LPStatus.UNBOUNDED; return; }

        pivot(r, e);
        status = LPStatus.RUNNING;
    }

    /** Entering: smallest column index with negative reduced cost (Bland-style order). */
    private int chooseEnteringColumn() {
        for (int c = 1; c <= n; c++) if (isNeg(T[0][c])) return c;
        return -1;
    }

    /**
     * Leaving row:
     *  - BLAND: ordinary min-ratio, tie by smallest row index.
     *  - LEXICOGRAPHIC: pick argmin of row/col-scaled vectors lexicographically:
     *      compare (T[r][0]/a_re, T[r][1]/a_re, ..., T[r][n]/a_re)
     *    across r with a_re > 0. This is the classic anti-cycling pivot lrs uses.
     */
    private int chooseLeavingRow(int enterCol) {
        if (pivotRule == PivotRule.BLAND) return leavingMinRatioBland(enterCol);
        return leavingLexicographic(enterCol);
    }

    private int leavingMinRatioBland(int e) {
        Fraction best = null;
        int arg = -1;
        for (int r = 1; r <= m; r++) {
            Fraction a = T[r][e];
            if (isPos(a)) {
                Fraction ratio = T[r][0].divide(a);
                if (best == null || ratio.compareTo(best) < 0) { best = ratio; arg = r; }
            }
        }
        return arg;
    }

    private int leavingLexicographic(int e) {
        int bestRow = -1;
        for (int r = 1; r <= m; r++) {
            if (!isPos(T[r][e])) continue;
            if (bestRow == -1) { bestRow = r; continue; }
            if (lexLess(r, bestRow, e)) bestRow = r;
        }
        return bestRow;
    }

    /**
     * Compare rows r vs s lexicographically after dividing each row by its pivot column (e):
     *   compare T[r][0]/T[r][e] vs T[s][0]/T[s][e], then T[r][1]/T[r][e] vs T[s][1]/T[s][e], ...
     * We avoid temporary fractions via cross-multiplication.
     */
    private boolean lexLess(int r, int s, int e) {
        Fraction re = T[r][e]; // > 0 guaranteed by caller
        Fraction se = T[s][e]; // > 0
        for (int c = 0; c <= n; c++) {
            int cmp = compareQuotients(T[r][c], re, T[s][c], se);
            if (cmp < 0) return true;
            if (cmp > 0) return false;
        }
        // Perfect tie (identical scaled rows) — fall back to smaller row index for determinism
        return r < s;
    }

    /** compare a/x ? b/y without constructing a/x and b/y: returns -1/0/1 */
    private static int compareQuotients(Fraction a, Fraction x, Fraction b, Fraction y) {
        // a/x ? b/y   <=>   a*y ? b*x
        return a.multiply(y).compareTo(b.multiply(x));
    }

    /** Pivot at (leaveRow, enterCol), Gauss–Jordan elimination and basis/cobasis swap. */
    public void pivot(int leaveRow, int enterCol) {
        Fraction piv = T[leaveRow][enterCol];
        if (isZero(piv)) throw new IllegalArgumentException("Pivot on zero element");

        // Normalize pivot row
        for (int c = 0; c <= n; c++) T[leaveRow][c] = T[leaveRow][c].divide(piv);

        // Eliminate column from other rows
        for (int r = 0; r <= m; r++) {
            if (r == leaveRow) continue;
            Fraction factor = T[r][enterCol];
            if (isZero(factor)) continue;
            for (int c = 0; c <= n; c++) {
                if (c == enterCol) {
                    T[r][c] = Fraction.ZERO;
                } else {
                    T[r][c] = T[r][c].subtract(factor.multiply(T[leaveRow][c]));
                }
            }
        }

        // Swap labels
        int enteringVar = cobasis[enterCol];
        int leavingVar  = basis[leaveRow];
        basis[leaveRow]   = enteringVar;
        cobasis[enterCol] = leavingVar;
    }

    /** Run until OPTIMAL/UNBOUNDED or iteration cap. */
    public LPStatus solve(int maxIters) {
        int it = 0;
        while (status == LPStatus.RUNNING && it < maxIters) { step(); it++; }
        return status;
    }
    enum PivotRule { BLAND, LEXICOGRAPHIC }
    enum LPStatus { RUNNING, OPTIMAL, UNBOUNDED }

    // --- helpers for reverse search (public on purpose)
    public int leavingRowFor(int enterCol) {
        // Uses the dictionary's active pivot rule (LEXICOGRAPHIC by default)
        // Returns -1 if no feasible leaving row (i.e., pivot would be infeasible/unbounded).
        return chooseLeavingRow(enterCol);
    }

    public boolean canPivot(int enterCol) {
        return leavingRowFor(enterCol) != -1;
    }

}
