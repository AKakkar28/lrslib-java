package com.vertexenumeration;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Locale;

/**
 * Represents an H‑ or V‑representation of a convex polyhedron.  This
 * class stores the header information (the type and field) and the
 * coefficient matrix.  It also provides static methods to read from
 * and write to the plain text file format used by lrslib.
 */
public class Polyhedron {
    /** Enumeration of representation types. */
    public enum Type { H, V }

    private final Type type;
    private final int rowCount;
    private final int colCount;
    private final Matrix matrix;
    private final boolean integerData;

    Polyhedron(Type type, int m, int n, boolean integerData, Matrix matrix) {
        this.type = type;
        this.rowCount = m;
        this.colCount = n;
        this.integerData = integerData;
        this.matrix = matrix;
    }

    public Type getType() { return type; }
    public int getRowCount() { return rowCount; }
    public int getColCount() { return colCount; }
    public boolean isIntegerData() { return integerData; }
    public Matrix getMatrix() { return matrix; }

    /**
     * Reads a polyhedron from an input file in lrslib format.  The
     * method consumes the entire file.  Example of accepted input:
     *
     * <pre>
     * H-representation
     * begin
     * 4 3 integer
     * 1 1 0
     * 1 0 1
     * 1 -1 0
     * 1 0 -1
     * end
     * </pre>
     *
     * @param filename path to the input file
     * @return a {@code Polyhedron} instance
     */
    public static Polyhedron readFromFile(String filename) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line;
            // Skip empty lines and comments
            do {
                line = br.readLine();
                if (line == null) {
                    throw new IOException("Unexpected end of file");
                }
                line = line.trim();
            } while (line.isEmpty());
            // Determine type (H- or V-representation)
            Type type;
            if (line.toLowerCase(Locale.ROOT).startsWith("h-representation")) {
                type = Type.H;
            } else if (line.toLowerCase(Locale.ROOT).startsWith("v-representation")) {
                type = Type.V;
            } else {
                throw new IOException("Expected H-representation or V-representation, got: " + line);
            }
            // Read 'begin'
            do {
                line = br.readLine();
                if (line == null) {
                    throw new IOException("Unexpected end of file");
                }
                line = line.trim();
            } while (line.isEmpty());
            if (!line.equalsIgnoreCase("begin")) {
                throw new IOException("Expected 'begin' after representation type, got: " + line);
            }
            // Read dimensions and integer/rational flag
            String[] header;
            do {
                line = br.readLine();
                if (line == null) {
                    throw new IOException("Unexpected end of file");
                }
                line = line.trim();
            } while (line.isEmpty());
            header = line.split("\\s+");
            if (header.length < 3) {
                throw new IOException("Expected 'm n integer|rational', got: " + line);
            }
            int m = Integer.parseInt(header[0]);
            int n = Integer.parseInt(header[1]);
            boolean integerData = header[2].equalsIgnoreCase("integer");
            Matrix matrix = Matrix.read(br, m, n);
            // Read 'end'
            do {
                line = br.readLine();
                if (line == null) {
                    throw new IOException("Unexpected end of file");
                }
                line = line.trim();
            } while (line.isEmpty());
            if (!line.equalsIgnoreCase("end")) {
                throw new IOException("Expected 'end', got: " + line);
            }
            return new Polyhedron(type, m, n, integerData, matrix);
        }
    }

    /**
     * Writes this polyhedron to the given {@link PrintWriter} in lrslib
     * format.  The caller is responsible for closing the writer.
     */
    public void write(PrintWriter out) throws IOException {
        if (type == Type.V) {
            // lrslib-style header for V-representation:
            out.println("V-representation");
            out.println("begin");
            // For V: lrslib prints "***** <cols> rational" (not m n).
            out.printf("***** %d rational%n", colCount);

            // Dump rows
            Matrix M = this.matrix;
            for (int i = 0; i < rowCount; i++) {
                for (int j = 0; j < colCount; j++) {
                    out.print(" " + M.get(i, j));
                }
                out.println();
            }
            out.println("end");
            return;
        }

        // H-representation stays as before
        out.println("H-representation");
        out.println("begin");
        out.printf("%d %d %s%n", rowCount, colCount, integerData ? "integer" : "rational");
        matrix.write(out);
        out.println("end");
    }

}