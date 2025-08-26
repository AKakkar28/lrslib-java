package com.vertexenumeration;

import java.util.Arrays;

/**
 * Build a canonical simplex dictionary from an H-representation:
 *
 *   H is m x (d+1): each row is [a0, a1, ..., a_d] meaning  a0 + a^T x >= 0.
 *
 * We construct a tableau with:
 *   - Row 0: objective (all zeros; we’re not solving an LP here)
 *   - Rows 1..m: slack-basis dictionary rows
 *       s_i = a0_i + sum_j a_{ij} x_j
 *
 * Columns (non-basic variables) are exactly x_1..x_d (so n = d).
 * Basis labels are the slacks s_1..s_m (we label them as d+1 .. d+m).
 */
final class HToDictionary {

    private HToDictionary() {}

    /**
     * @param H m x (d+1) matrix where H[i][0] = a0_i, H[i][j] = a_{i,j} for j>=1.
     * @return Dictionary with slacks basic and x’s non-basic.
     */
    static Dictionary fromH(Fraction[][] H) {
        if (H == null || H.length == 0) {
            throw new IllegalArgumentException("H must be non-empty");
        }
        final int m = H.length;
        final int ncols = H[0].length;
        if (ncols < 2) {
            throw new IllegalArgumentException("H must have at least a0 and one variable column");
        }
        final int d = ncols - 1; // number of x-variables

        // Sanity: rectangular
        for (int i = 0; i < m; i++) {
            if (H[i].length != ncols) {
                throw new IllegalArgumentException("All rows of H must have same length");
            }
            for (int j = 0; j < ncols; j++) {
                if (H[i][j] == null) throw new IllegalArgumentException("Null Fraction in H at ("+i+","+j+")");
            }
        }

        // Build tableau T of shape (m+1) x (d+1).
        // T[0][*] is objective row: everything ZERO.
        // T[i][0] = a0_i (RHS)
        // T[i][j] = a_{i,j}  for j=1..d
        Fraction ZERO = H[0][0].subtract(H[0][0]);
        Fraction[][] T = new Fraction[m + 1][d + 1];

        // Objective row (zeros)
        Arrays.fill(T[0] , ZERO);

        // Constraint rows
        for (int i = 1; i <= m; i++) {
            T[i][0] = H[i - 1][0];            // RHS = a0_i
            for (int j = 1; j <= d; j++) {
                T[i][j] = H[i - 1][j];        // coeff of x_j
            }
        }

        // Basis = slacks s_1..s_m labeled as d+1..d+m
        int[] basis = new int[m + 1];
        basis[0] = 0; // unused
        for (int i = 1; i <= m; i++) basis[i] = d + i;

        // Cobasis = structural variables x_1..x_d labeled 1..d
        int[] cobasis = new int[d + 1];
        cobasis[0] = 0; // unused
        for (int j = 1; j <= d; j++) cobasis[j] = j;

        return Dictionary.of(T, basis, cobasis);
    }
}
