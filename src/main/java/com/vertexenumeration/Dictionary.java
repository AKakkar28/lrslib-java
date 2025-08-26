package com.vertexenumeration;

import java.util.Arrays;
import java.util.Objects;

/**
 * Dictionary (simplex tableau) for reverse search / LP steps.
 *
 * Layout:
 *   Row 0   : objective (z - c^T x = 0), RHS in col 0
 *   Rows 1..m : constraints (b - A x = 0), RHS in col 0
 *   Col 0   : RHS
 *   Cols 1..n: non-basic variable columns
 *
 * basis[r]   = variable id basic in row r (1..m)
 * cobasis[c] = variable id non-basic in column c (1..n)
 */
final class Dictionary {

    private final int m;            // constraints
    private final int n;            // non-basic columns
    private final Fraction[][] T;   // (m+1) x (n+1) tableau; T[row][col]
    private final int[] basis;      // size m+1, index 0 unused
    private final int[] cobasis;    // size n+1, index 0 unused
    private PivotRule pivotRule = PivotRule.BLAND;

    private LPStatus status = LPStatus.RUNNING;

    // Local ZERO derived from the incoming tableau; avoids relying on Fraction.ZERO
    private final Fraction ZERO;

    Dictionary(int m, int n, Fraction[][] tableau, int[] basis, int[] cobasis) {
        this.m = m;
        this.n = n;
        this.T = deepCopy(tableau);
        this.basis = Arrays.copyOf(basis, basis.length);
        this.cobasis = Arrays.copyOf(cobasis, cobasis.length);
        sanity();

        // Derive a ZERO value without assuming a static constant on Fraction
        Fraction probe = this.T[0][0];
        this.ZERO = probe.subtract(probe);
    }

    // Convenience builder
    public static Dictionary of(Fraction[][] tableau, int[] basis, int[] cobasis) {
        int m = tableau.length - 1;
        int n = tableau[0].length - 1;
        return new Dictionary(m, n, tableau, basis, cobasis);
    }

    private void sanity() {
        if (T.length != m + 1) throw new IllegalArgumentException("tableau row count mismatch");
        for (Fraction[] row : T) {
            if (row.length != n + 1) throw new IllegalArgumentException("tableau col count mismatch");
            for (int j = 0; j <= n; j++) {
                if (row[j] == null) throw new IllegalArgumentException("null Fraction in tableau");
            }
        }
        if (basis.length != m + 1) throw new IllegalArgumentException("basis length mismatch");
        if (cobasis.length != n + 1) throw new IllegalArgumentException("cobasis length mismatch");
    }

    private static Fraction[][] deepCopy(Fraction[][] a) {
        Fraction[][] c = new Fraction[a.length][];
        for (int i = 0; i < a.length; i++) c[i] = Arrays.copyOf(a[i], a[i].length);
        return c;
    }

    // --- accessors (useful for your ReverseSearchEnumerator wiring)
    public int m() { return m; }
    public int n() { return n; }
    public LPStatus status() { return status; }
    public int[] basis() { return Arrays.copyOf(basis, basis.length); }
    public int[] cobasis() { return Arrays.copyOf(cobasis, cobasis.length); }
    public Fraction[][] tableau() { return deepCopy(T); }

    public void setPivotRule(PivotRule rule) { this.pivotRule = Objects.requireNonNull(rule); }

    // --- small helpers that don’t rely on Fraction methods
    private boolean isNeg(Fraction x) { return x.compareTo(ZERO) < 0; }
    private boolean isPos(Fraction x) { return x.compareTo(ZERO) > 0; }
    private boolean isZero(Fraction x){ return x.compareTo(ZERO) == 0; }

    /** Feasibility for standard form: all RHS (rows 1..m, col 0) >= 0. */
    public boolean isFeasible() {
        for (int r = 1; r <= m; r++) if (isNeg(T[r][0])) return false;
        return true;
    }

    /** Any negative reduced cost in row 0, cols 1..n? */
    private boolean anyNegativeReducedCost() {
        for (int c = 1; c <= n; c++) if (isNeg(T[0][c])) return true;
        return false;
    }

    /**
     * One simplex step:
     *  - choose entering column (Bland)
     *  - ratio test for leaving row
     *  - pivot
     * Sets status to OPTIMAL / UNBOUNDED when appropriate.
     */
    public void step() {
        if (!anyNegativeReducedCost()) { status = LPStatus.OPTIMAL; return; }

        int enterCol = chooseEnteringColumn();
        if (enterCol == -1) { status = LPStatus.OPTIMAL; return; }

        int leaveRow = chooseLeavingRow(enterCol);
        if (leaveRow == -1) { status = LPStatus.UNBOUNDED; return; }

        pivot(leaveRow, enterCol);
        status = LPStatus.RUNNING;
    }

    /** Bland’s rule: smallest column index with negative reduced cost. */
    private int chooseEnteringColumn() {
        if (pivotRule == PivotRule.BLAND) {
            for (int c = 1; c <= n; c++) if (isNeg(T[0][c])) return c;
            return -1;
        }
        // Future: DANTZIG, STEEPEST_EDGE, etc.
        for (int c = 1; c <= n; c++) if (isNeg(T[0][c])) return c;
        return -1;
    }

    /**
     * Min ratio test: argmin { RHS / a_rc | a_rc > 0 } across rows r=1..m.
     * @return leaving row index or -1 if a_rc <= 0 for all r (unbounded).
     */
    private int chooseLeavingRow(int enterCol) {
        Fraction best = null;
        int arg = -1;
        for (int r = 1; r <= m; r++) {
            Fraction a = T[r][enterCol];
            if (isPos(a)) {
                Fraction ratio = T[r][0].divide(a);
                if (best == null || ratio.compareTo(best) < 0) {
                    best = ratio;
                    arg = r;
                }
            }
        }
        return arg;
    }

    /**
     * Pivot at (leaveRow, enterCol). Gauss–Jordan update and label swap.
     */
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
                    T[r][c] = ZERO; // becomes exactly zero
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

    /**
     * Run simplex until OPTIMAL/UNBOUNDED or iteration cap.
     * Returns final status (RUNNING can mean degeneracy protection hit maxIters).
     */
    public LPStatus solve(int maxIters) {
        int it = 0;
        while (status == LPStatus.RUNNING && it < maxIters) {
            step();
            it++;
        }
        return status;
    }
}
