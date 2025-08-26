package com.vertexenumeration;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Represents an H- or V-representation of a convex polyhedron.
 * Parses (and writes) an lrslib-compatible format. Output uses the lrs
 * style line: "***** <ncols> rational" for both H and V.
 */
public class Polyhedron {
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

    public static Polyhedron readFromFile(String filename) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line;

            // ---- scan header (allow: name line, options, then representation/begin) ----
            Type type = Type.H;               // default if omitted (lrs behavior)
            boolean sawBegin = false;

            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.isEmpty()) continue;
                if (line.startsWith("*") || line.startsWith("#")) continue; // comment

                String low = line.toLowerCase(Locale.ROOT);
                if (low.startsWith("h-representation")) { type = Type.H; continue; }
                if (low.startsWith("v-representation")) { type = Type.V; continue; }

                if (low.equals("begin")) { sawBegin = true; break; }

                // otherwise treat as a free "name"/option line and ignore
            }
            if (!sawBegin) throw new IOException("No 'begin' line found");

            // ---- read header line after 'begin' ----
            String[] header;
            do {
                line = br.readLine();
                if (line == null) throw new IOException("Unexpected end of file");
                line = line.trim();
            } while (line.isEmpty());
            header = line.split("\\s+");
            if (header.length < 3)
                throw new IOException("Expected 'm n integer|rational' or '***** n integer|rational', got: " + line);

            final boolean starHeader = header[0].equals("*****");
            final int n = Integer.parseInt(header[1]);
            final boolean integerData = header[2].equalsIgnoreCase("integer");

            Matrix matrix;
            int m;

            if (starHeader) {
                // Read rows until 'end'
                List<String> rows = new ArrayList<>();
                while ((line = br.readLine()) != null) {
                    String t = line.trim();
                    if (t.isEmpty()) continue;
                    if (t.equalsIgnoreCase("end")) break;
                    rows.add(t);
                }
                m = rows.size();

                // Re-read the collected rows through a StringReader so Matrix.read can parse them
                StringBuilder sb = new StringBuilder();
                for (String r : rows) sb.append(r).append('\n');
                try (BufferedReader tmp = new BufferedReader(new StringReader(sb.toString()))) {
                    matrix = Matrix.read(tmp, m, n);
                }
                // 'end' already consumed above
            } else {
                // Numeric row count
                if (!header[0].matches("[-+]?\\d+"))
                    throw new IOException("Row count must be numeric or '*****', got: " + header[0]);
                m = Integer.parseInt(header[0]);
                matrix = Matrix.read(br, m, n);

                // Expect 'end'
                do {
                    line = br.readLine();
                    if (line == null) throw new IOException("Unexpected end of file");
                    line = line.trim();
                } while (line.isEmpty());
                if (!line.equalsIgnoreCase("end"))
                    throw new IOException("Expected 'end', got: " + line);
            }

            return new Polyhedron(type, m, n, integerData, matrix);
        }
    }

    /** Write using lrs-style starred header for both H and V. */
    public void write(PrintWriter out) throws IOException {
        out.println((type == Type.H ? "H" : "V") + "-representation");
        out.println("begin");
        out.printf("***** %d %s%n", colCount, "rational"); // lrs-style
        matrix.write(out);
        out.println("end");
    }
}
