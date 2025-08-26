package com.vertexenumeration;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Simple command‑line front‑end for the vertex enumeration program.  It
 * reads an input file in lrslib format, converts it to the dual
 * representation, and writes the result to standard output.
 */
public class Main {
    public static void main(String[] args) {
        if (args.length != 1) {
            System.err.println("Usage: java -jar vertex-enumeration-java.jar <input-file>");
            System.exit(1);
        }
        String filename = args[0];
        try {
            Polyhedron input = Polyhedron.readFromFile(filename);
            VertexEnumerator enumerator = new VertexEnumerator();
            Polyhedron output = enumerator.enumerate(input);
            PrintWriter pw = new PrintWriter(System.out);
            output.write(pw);
            pw.flush();
        } catch (FileNotFoundException e) {
            System.err.println("File not found: " + filename);
            System.exit(1);
        } catch (IOException e) {
            System.err.println("I/O error: " + e.getMessage());
            System.exit(1);
        }
    }
}