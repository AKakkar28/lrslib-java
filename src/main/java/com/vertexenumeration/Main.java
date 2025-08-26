package com.vertexenumeration;

import java.io.*;

public class Main {
    public static void main(String[] args) {
        if (args.length != 1) {
            System.err.println("Usage: java -cp target/<jar> com.vertexenumeration.Main <input-file>");
            System.exit(1);
        }
        String filename = args[0];
        long t0 = System.nanoTime();

        try {
            Polyhedron input = Polyhedron.readFromFile(filename);
            VertexEnumerator enumerator = new VertexEnumerator();
            Polyhedron output = enumerator.enumerate(input);

            // write output
            PrintWriter pw = new PrintWriter(System.out);
            output.write(pw);
            pw.flush();

            // footer (*Totals...) like lrs
            EnumStats st = enumerator.getLastStats();
            long t1 = System.nanoTime();
            if (st != null) {
                System.out.println(st.toString());
            }
            // simple timing line (we won’t mimic lrs’s exact format yet)
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
