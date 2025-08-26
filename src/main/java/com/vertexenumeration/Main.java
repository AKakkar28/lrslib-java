package com.vertexenumeration;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Paths;

public class Main {
    public static void main(String[] args) {
        if (args.length != 1) {
            System.err.println("Usage: java -cp target/<jar or classes> com.vertexenumeration.Main <input-file>");
            System.exit(1);
        }
        String filename = args[0];
        long t0 = System.nanoTime();

        try {
            Polyhedron input = Polyhedron.readFromFile(filename);
            VertexEnumerator enumerator = new VertexEnumerator();
            Polyhedron output = enumerator.enumerate(input);

            // lrs-style "name" line (base filename w/o extension)
            String base = Paths.get(filename).getFileName().toString();
            int dot = base.lastIndexOf('.');
            if (dot > 0) base = base.substring(0, dot);
            System.out.println(base);

            // body
            PrintWriter pw = new PrintWriter(System.out);
            output.write(pw);
            pw.flush();

            // footer
            EnumStats st = enumerator.getLastStats();
            long t1 = System.nanoTime();
            if (st != null) {
                if (output.getType() == Polyhedron.Type.V) {
                    // H -> V case: your existing Stats string
                    System.out.println(st.toString());
                } else {
                    // V -> H case: report facets like lrs
                    System.out.printf("*Totals: facets=%d bases=%d%n",
                            output.getRowCount(), st.bases);
                }
            }
            double secs = (t1 - t0) / 1_000_000_000.0;
            System.out.printf("*Time=%.3fs%n", secs);

        } catch (FileNotFoundException e) {
            System.err.println("File not found: " + filename);
            System.exit(1);
        } catch (IOException e) {
            System.err.println("I/O error: " + e.getMessage());
            System.exit(1);
        }
    }
}
