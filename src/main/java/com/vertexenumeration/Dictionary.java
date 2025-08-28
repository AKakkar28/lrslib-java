package com.vertexenumeration;

import java.util.Arrays;
import java.util.Objects;

/**
 * Dictionary structure with full parity to lrslib's lrs_dic_struct.
 * Maintains tableau, basis/cobasis, row/col permutations, and pivot state.
 */
final class Dictionary implements LrsDic {

    private final int m, n;             // constraints, variables
    private final Fraction[][] T;       // tableau (m+1) x (n+1)

    // Basis and cobasis arrays
    private final int[] basis;          // size m+1, basic vars
    private final int[] cobasis;        // size n+1, nonbasic vars

    // Row/Col permutations (for exact C parity)
    private final int[] Row;            // Row[r] = original row index at position r
    private final int[] Col;            // Col[c] = original col index at position c

    // Explicit copies of basis/nonbasis (like C’s B[], C[])
    private final int[] B;              // basic var indices
    private final int[] C;              // nonbasic var indices

    // State fields
    private int depth;                  // recursion depth
    private int pivotRow, pivotCol;     // last pivot position
    private int lastDecisionVar;        // last decision var index
    private int objCol;                 // objective column index

    private PivotRule pivotRule = PivotRule.LEXICOGRAPHIC;
    private LPStatus status = LPStatus.RUNNING;

    Dictionary(int m, int n, Fraction[][] tableau, int[] basis, int[] cobasis) {
        this.m = m;
        this.n = n;
        this.T = deepCopy(tableau);

        this.basis = Arrays.copyOf(basis, basis.length);
        this.cobasis = Arrays.copyOf(cobasis, cobasis.length);

        this.Row = new int[m + 1];
        this.Col = new int[n + 1];
        for (int i = 0; i <= m; i++) Row[i] = i;
        for (int j = 0; j <= n; j++) Col[j] = j;

        this.B = Arrays.copyOf(basis, basis.length);
        this.C = Arrays.copyOf(cobasis, cobasis.length);

        this.depth = 0;
        this.pivotRow = -1;
        this.pivotCol = -1;
        this.lastDecisionVar = -1;
        this.objCol = 0; // default objective column index

        sanity();
    }

    // --- Accessors
    @Override
    public int m() { return m; }
    @Override
    public int n() { return n; }
    @Override
    public Fraction[][] tableau() { return deepCopy(T); }
    @Override
    public int[] basis() { return Arrays.copyOf(basis, basis.length); }
    @Override
    public int[] cobasis() { return Arrays.copyOf(cobasis, cobasis.length); }
    @Override
    public LPStatus status() { return status; }

    public int[] rowPerm() { return Arrays.copyOf(Row, Row.length); }
    public int[] colPerm() { return Arrays.copyOf(Col, Col.length); }
    public int[] B() { return Arrays.copyOf(B, B.length); }
    public int[] C() { return Arrays.copyOf(C, C.length); }
    public int depth() { return depth; }
    public int lastDecisionVar() { return lastDecisionVar; }
    public int objCol() { return objCol; }

    @Override
    public boolean isFeasible() {
        for (int r = 1; r <= m; r++) if (T[r][0].compareTo(Fraction.ZERO) < 0) return false;
        return true;
    }

    @Override
    public void setPivotRule(PivotRule rule) { this.pivotRule = Objects.requireNonNull(rule); }

    @Override
    public int leavingRowFor(int enterCol) { return chooseLeavingRow(enterCol); }

    @Override
    public boolean canPivot(int enterCol) { return leavingRowFor(enterCol) != -1; }

    @Override
    public LPStatus solve(int maxIters) {
        int it = 0;
        while (status == LPStatus.RUNNING && it < maxIters) {
            step();
            it++;
        }
        return status;
    }

    // --- Simplex iteration logic
    private void step() {
        if (!anyNegativeReducedCost()) { status = LPStatus.OPTIMAL; return; }
        int e = chooseEnteringColumn();
        if (e == -1) { status = LPStatus.OPTIMAL; return; }
        int r = chooseLeavingRow(e);
        if (r == -1) { status = LPStatus.UNBOUNDED; return; }
        pivot(r, e);
        status = LPStatus.RUNNING;
    }

    private boolean anyNegativeReducedCost() {
        for (int c = 1; c <= n; c++) if (T[0][c].compareTo(Fraction.ZERO) < 0) return true;
        return false;
    }

    private int chooseEnteringColumn() {
        for (int c = 1; c <= n; c++) if (T[0][c].compareTo(Fraction.ZERO) < 0) return c;
        return -1;
    }

    private int chooseLeavingRow(int enterCol) {
        if (pivotRule == PivotRule.BLAND) return leavingMinRatioBland(enterCol);
        return leavingLexicographic(enterCol);
    }

    private int leavingMinRatioBland(int e) {
        Fraction best = null;
        int arg = -1;
        for (int r = 1; r <= m; r++) {
            Fraction a = T[r][e];
            if (a.compareTo(Fraction.ZERO) > 0) {
                Fraction ratio = T[r][0].divide(a);
                if (best == null || ratio.compareTo(best) < 0) { best = ratio; arg = r; }
            }
        }
        return arg;
    }

    private int leavingLexicographic(int e) {
        int bestRow = -1;
        for (int r = 1; r <= m; r++) {
            Fraction a = T[r][e];
            if (a.compareTo(Fraction.ZERO) <= 0) continue;
            if (bestRow == -1) { bestRow = r; continue; }
            if (lexLess(r, bestRow, e)) bestRow = r;
        }
        return bestRow;
    }

    private boolean lexLess(int r, int s, int e) {
        Fraction re = T[r][e], se = T[s][e];
        for (int c = 0; c <= n; c++) {
            int cmp = compareQuotients(T[r][c], re, T[s][c], se);
            if (cmp < 0) return true;
            if (cmp > 0) return false;
        }
        return r < s;
    }

    private static int compareQuotients(Fraction a, Fraction x, Fraction b, Fraction y) {
        return a.multiply(y).compareTo(b.multiply(x));
    }

    @Override
    public void pivot(int leaveRow, int enterCol) {
        Fraction piv = T[leaveRow][enterCol];
        if (piv.compareTo(Fraction.ZERO) == 0) throw new IllegalArgumentException("Pivot on zero element");

        for (int c = 0; c <= n; c++) T[leaveRow][c] = T[leaveRow][c].divide(piv);

        for (int r = 0; r <= m; r++) {
            if (r == leaveRow) continue;
            Fraction factor = T[r][enterCol];
            if (factor.isZero()) continue;
            for (int c = 0; c <= n; c++) {
                if (c == enterCol) T[r][c] = Fraction.ZERO;
                else T[r][c] = T[r][c].subtract(factor.multiply(T[leaveRow][c]));
            }
        }

        int enteringVar = cobasis[enterCol];
        int leavingVar  = basis[leaveRow];
        basis[leaveRow]   = enteringVar;
        cobasis[enterCol] = leavingVar;

        B[leaveRow] = enteringVar;
        C[enterCol] = leavingVar;
        Row[leaveRow] = enteringVar;
        Col[enterCol] = leavingVar;

        this.pivotRow = leaveRow;
        this.pivotCol = enterCol;
        this.lastDecisionVar = enteringVar;
        this.depth++;
    }

    private static Fraction[][] deepCopy(Fraction[][] a) {
        Fraction[][] c = new Fraction[a.length][];
        for (int i = 0; i < a.length; i++) c[i] = Arrays.copyOf(a[i], a[i].length);
        return c;
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

    enum PivotRule { BLAND, LEXICOGRAPHIC }
    enum LPStatus { RUNNING, OPTIMAL, UNBOUNDED }
}