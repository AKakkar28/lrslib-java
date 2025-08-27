package com.vertexenumeration;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Paths;

public class Main {

    private static void usage() {
        System.err.println(
                "Usage: lrs-java [options] <input-file>\n" +
                        "Modes:\n" +
                        "  -v            H->V (vertex enumeration)  [default]\n" +
                        "  -h            V->H (facet enumeration)\n" +
                        "Options:\n" +
                        "  -redund       remove redundant inequalities (H only)\n" +
                        "  -minrep       hidden linearities + minimal H (H only)\n" +
                        "  -eliminate c1,c2,...   Fourier–Motzkin eliminate columns (H only)\n" +
                        "  -project   c1,c2,...   keep only these columns (H only)\n" +
                        "  -linearity r1,r2,...   mark rows as equalities (H only)\n" +
                        "  -printcobasis           print cobasis indices with output (todo)\n" +
                        "  -seed N       deterministic tie-breaking seed\n" +
                        "  -maxdepth N    cap reverse-search depth (0 = no cap)\n" +
                        "  -threads N     root parallelism (implementation-defined)\n" +
                        "  -integer       declare input as integer data (affects stats only)\n"
        );
    }

    public static void main(String[] args) {
        OptionsParser.Parsed parsed;
        try {
            parsed = OptionsParser.parse(args);
        } catch (Exception e) {
            usage();
            System.err.println("Argument error: " + e.getMessage());
            System.exit(2);
            return;
        }

        LrsDat dat = parsed.dat;
        String filename = parsed.inputPath;

        long t0 = System.nanoTime();

        try {
            // read
            Polyhedron input = Polyhedron.readFromFile(filename);

            // apply Step-2 transform hooks (you’ll implement these next)
            Polyhedron working = Transforms.applyAll(input, dat);

            // run
            Polyhedron output;
            EnumStats st = null;

            if (dat.mode == LrsDat.Mode.VE) {
                // H -> V
                VertexEnumerator enumerator = new VertexEnumerator();
                output = enumerator.enumerate(working);
                st = enumerator.getLastStats();
                System.err.println("DEBUG: stats vertices=" + st.vertices);

            } else {
                // V -> H
                output = FacetEnumerator.fromV(working);
                // we don't compute detailed stats for V->H yet
                st = new EnumStats();
                st.bases = output.getRowCount(); // placeholder like before
            }

            // lrs-style "name" line (base filename without extension)
            String base = Paths.get(filename).getFileName().toString();
            int dot = base.lastIndexOf('.');
            if (dot > 0) base = base.substring(0, dot);
            System.out.println(base);

            // body
            PrintWriter pw = new PrintWriter(System.out, true);
            output.write(pw);
            pw.flush();

            long t1 = System.nanoTime();
            if (st != null) {
                if (dat.mode == LrsDat.Mode.VE) {
                    System.out.printf("*Totals: vertices=%d rays=%d bases=%d integer_vertices=%d%n",
                            st.vertices, st.rays, st.bases, st.integerVertices);
                    System.out.printf("*max_vertex_depth=%d%n", st.maxDepth);
                } else {
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
