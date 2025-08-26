package com.vertexenumeration;

import java.util.ArrayList;
import java.util.List;

/**
 * Represents a simplex dictionary (tableau) used in the reverse search
 * algorithm for vertex enumeration.  Each dictionary corresponds to a
 * particular basis/cobasis and stores the coefficients of the linear
 * constraints and objective.  This class only defines the data
 * structures; the pivot operations and lexicographic rules must be
 * implemented separately.
 */
public class Dictionary {
    /** Number of constraints (rows). */
    public final int m;
    /** Number of variables (columns). */
    public final int n;

    /** The tableau coefficients.  An (m×n) matrix of Fractions. */
    public final Fraction[][] tableau;

    /** Indices of basic variables.  Length m. */
    public final int[] basis;

    /** Indices of non-basic variables.  Length n. */
    public final int[] cobasis;

    /** Current solution vector for the basic variables.  Length m. */
    public final Fraction[] solution;

    /**
     * Constructs a blank dictionary of the given dimensions.  All
     * fractional entries are initialised to zero.  Basis and cobasis
     * arrays are initialised with default indices (0..m-1 for the
     * basis, m..m+n-1 for the cobasis).
     */
    public Dictionary(int m, int n) {
        this.m = m;
        this.n = n;
        this.tableau = new Fraction[m][n];
        Fraction zero = new Fraction(java.math.BigInteger.ZERO, java.math.BigInteger.ONE);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                tableau[i][j] = zero;
            }
        }
        this.basis = new int[m];
        this.cobasis = new int[n];
        this.solution = new Fraction[m];
        for (int i = 0; i < m; i++) {
            basis[i] = i;
            solution[i] = zero;
        }
        for (int j = 0; j < n; j++) {
            cobasis[j] = m + j;
        }
    }

    /**
     * Returns a copy of this dictionary.  Only the arrays are cloned; the
     * Fraction objects themselves are immutable and thus shared.
     */
    public Dictionary copy() {
        Dictionary d = new Dictionary(m, n);
        for (int i = 0; i < m; i++) {
            d.basis[i] = this.basis[i];
            d.solution[i] = this.solution[i];
            for (int j = 0; j < n; j++) {
                d.tableau[i][j] = this.tableau[i][j];
            }
        }
        for (int j = 0; j < n; j++) {
            d.cobasis[j] = this.cobasis[j];
        }
        return d;
    }

    /**
     * Performs a pivot operation on the dictionary.  Given a pivot row
     * {@code r} and a pivot column {@code c}, this method returns a new
     * dictionary corresponding to the adjacent basis obtained by
     * entering variable c and leaving variable r.  **Note:** this is a
     * placeholder implementation that does not yet enforce the
     * lexicographic pivot rule or check feasibility.  Proper pivot
     * operations must be implemented when translating the algorithm.
     *
     * @param r the index of the leaving variable (row)
     * @param c the index of the entering variable (column)
     * @return a new {@code Dictionary} resulting from the pivot
     */
    public Dictionary pivot(int r, int c) {
        Dictionary next = this.copy();
        // Update basis and cobasis indices
        int leavingVar = next.basis[r];
        int enteringVar = next.cobasis[c];
        next.basis[r] = enteringVar;
        next.cobasis[c] = leavingVar;
        // Compute new tableau coefficients
        // Pivot element
        Fraction pivot = tableau[r][c];
        if (pivot.compareTo(new Fraction(java.math.BigInteger.ZERO, java.math.BigInteger.ONE)) == 0) {
            throw new ArithmeticException("Pivot element is zero");
        }
        // Update pivot row
        for (int j = 0; j < n; j++) {
            next.tableau[r][j] = tableau[r][j].divide(pivot);
        }
        next.solution[r] = solution[r].divide(pivot);
        // Update remaining rows
        for (int i = 0; i < m; i++) {
            if (i == r) continue;
            Fraction factor = tableau[i][c];
            for (int j = 0; j < n; j++) {
                Fraction val = tableau[i][j].subtract(factor.multiply(next.tableau[r][j]));
                next.tableau[i][j] = val;
            }
            Fraction solVal = solution[i].subtract(factor.multiply(next.solution[r]));
            next.solution[i] = solVal;
        }
        // Set pivot column in pivot row and pivot column in other rows
        for (int i = 0; i < m; i++) {
            next.tableau[i][c] = (i == r)
                    ? new Fraction(java.math.BigInteger.ONE, java.math.BigInteger.ONE)
                    : new Fraction(java.math.BigInteger.ZERO, java.math.BigInteger.ONE);
        }
        // Set pivot row solution entry to pivot row solution
        next.solution[r] = solution[r].divide(pivot);
        return next;
    }
}