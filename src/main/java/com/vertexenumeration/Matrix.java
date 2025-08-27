package com.vertexenumeration;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.io.Writer;
import java.math.BigInteger;

/**
 * A simple dense matrix of {@link Fraction} values.  Rows and columns are
 * indexed from zero.  This class is intentionally minimal; additional
 * operations can be added as needed for the vertex enumeration algorithm.
 */
public class Matrix {
    private final int rows;
    private final int cols;
    private final Fraction[][] data;

    public int rows() { return rows; }
    public int cols() { return cols; }
    /**
     * Constructs a {@code rows × cols} matrix with all entries set to
     * zero.
     */
    public Matrix(int rows, int cols) {
        if (rows < 0 || cols < 0) {
            throw new IllegalArgumentException("Negative dimensions");
        }
        this.rows = rows;
        this.cols = cols;
        this.data = new Fraction[rows][cols];
        Fraction zero = new Fraction(BigInteger.ZERO, BigInteger.ONE);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                data[i][j] = zero;
            }
        }
    }

    /** Returns the number of rows. */
    public int getRowCount() { return rows; }

    /** Returns the number of columns. */
    public int getColumnCount() { return cols; }

    /** Returns the entry at row {@code r}, column {@code c}. */
    public Fraction get(int r, int c) {
        return data[r][c];
    }

    /** Sets the entry at row {@code r}, column {@code c}. */
    public void set(int r, int c, Fraction value) {
        data[r][c] = value;
    }

    /**
     * Reads a matrix from an lrslib-style list of rows.  Each line should
     * contain {@code cols} whitespace‑separated integers or fractions.  A
     * fraction must be in the form {@code numerator/denominator}.  The
     * caller is responsible for reading the header lines (such as
     * “H‑representation begin m n integer”) before invoking this method.
     */
    public static Matrix read(BufferedReader reader, int m, int n) throws IOException {
        Matrix mat = new Matrix(m, n);
        for (int i = 0; i < m; i++) {
            String line;
            do {
                line = reader.readLine();
                if (line == null) {
                    throw new IOException("Unexpected end of file while reading matrix");
                }
                line = line.trim();
            } while (line.isEmpty());
            String[] tokens = line.split("\\s+");
            if (tokens.length != n) {
                throw new IOException("Expected " + n + " columns on line " + (i + 1));
            }
            for (int j = 0; j < n; j++) {
                mat.set(i, j, parseFraction(tokens[j]));
            }
        }
        return mat;
    }

    /**
     * Writes this matrix to a writer in lrslib form, one row per line.
     */
    public void write(Writer out) throws IOException {
        for (int i = 0; i < rows; i++) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < cols; j++) {
                if (j > 0) sb.append(' ');
                sb.append(data[i][j].toString());
            }
            out.write(sb.toString());
            out.write(System.lineSeparator());
        }
    }

    private static Fraction parseFraction(String token) {
        if (token.contains("/")) {
            String[] parts = token.split("/");
            BigInteger num = new BigInteger(parts[0]);
            BigInteger den = new BigInteger(parts[1]);
            return new Fraction(num, den);
        } else {
            return new Fraction(new BigInteger(token), BigInteger.ONE);
        }
    }
}