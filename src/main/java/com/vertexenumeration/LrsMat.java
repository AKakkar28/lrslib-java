package com.vertexenumeration;

import java.util.Objects;

public final class LrsMat {
    private final Matrix M;

    public LrsMat(Matrix backing) { this.M = Objects.requireNonNull(backing); }

    public int rows() { return M.rows(); }
    public int cols() { return M.cols(); }

    public Fraction get(int r, int c) { return M.get(r, c); }
    public void set(int r, int c, Fraction v) { M.set(r, c, v); }

    public RowView row(int r) { return new RowView(r); }
    public ColView col(int c) { return new ColView(c); }

    public LrsMat submatrix(int r0, int r1, int c0, int c1) {   // note: “submatrix” spelling
        int R = r1 - r0, C = c1 - c0;
        Matrix S = new Matrix(R, C);
        for (int i = 0; i < R; i++)
            for (int j = 0; j < C; j++)
                S.set(i, j, M.get(r0 + i, c0 + j));
        return new LrsMat(S);
    }

    public Fraction[][] toArray() {
        Fraction[][] A = new Fraction[rows()][cols()];
        for (int i = 0; i < rows(); i++)
            for (int j = 0; j < cols(); j++)
                A[i][j] = M.get(i, j);
        return A;
    }

    public final class RowView {
        private final int r;
        private RowView(int r) { this.r = r; }
        public int length() { return cols(); }
        public Fraction get(int j) { return M.get(r, j); }
        public void set(int j, Fraction v) { M.set(r, j, v); }
        public Fraction[] toArray() {
            Fraction[] a = new Fraction[cols()];
            for (int j = 0; j < cols(); j++) a[j] = M.get(r, j);
            return a;
        }
    }

    public final class ColView {
        private final int c;
        private ColView(int c) { this.c = c; }
        public int length() { return rows(); }
        public Fraction get(int i) { return M.get(i, c); }
        public void set(int i, Fraction v) { M.set(i, c, v); }
        public Fraction[] toArray() {
            Fraction[] a = new Fraction[rows()];
            for (int i = 0; i < rows(); i++) a[i] = M.get(i, c);
            return a;
        }
    }
}
