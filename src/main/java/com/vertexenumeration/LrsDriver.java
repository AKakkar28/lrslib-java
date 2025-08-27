package com.vertexenumeration;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.Duration;
import java.time.Instant;
import java.util.Arrays;

/**
 * Java driver that mirrors the structure of lrslib's lrsdriver.c / lrs.c.
 *
 * Responsibilities (aligned with the C driver):
 *  - Parse CLI (delegates to OptionsParser you already have)
 *  - Set up restart/state bookkeeping (RestartData)
 *  - Read the polyhedron input file
 *  - Apply requested transforms (delegates to Transforms)
 *  - Dispatch to the appropriate enumeration mode (H->V or V->H)
 *  - Print results, stats, and elapsed time
 *
 * Notes on arithmetic restarts:
 *  lrslib can restart with wider integer types (128-bit or GMP) when overflow is suspected.
 *  In Java we use BigInteger-backed Fraction, so arithmetic overflow restarts are unnecessary.
 *  We still keep the RestartData scaffolding for feature parity and future extensions
 *  (e.g., partial enumerations, checkpoints).
 */
public final class LrsDriver {

    /** Entry point to run one lrs job, similar to lrs_main(). */
    public int run(String[] args) {
        OptionsParser.Parsed parsed;
        try {
            parsed = OptionsParser.parse(args);
        } catch (Exception e) {
            MainErrorPrinter.usage();
            System.err.println("Argument error: " + e.getMessage());
            return 2; // match typical CLI convention
        }

        final LrsDat dat = parsed.dat;
        final String inputPath = parsed.inputPath;

        final RestartData R = new RestartData();
        R.resetDefaults();

        final Instant t0 = Instant.now();
        try {
            // 1) Read input polyhedron
            Polyhedron input = Polyhedron.readFromFile(inputPath);

            // 2) Apply requested transforms (projection, linearities, elimination, etc.)
            Polyhedron working = Transforms.applyAll(input, dat);

            // 3) Dispatch to algorithm
            Polyhedron output;
            EnumStats stats = null;

            switch (dat.mode) {
                case VE: { // H -> V (vertex/ray enumeration)
                    VertexEnumerator enumerator = new VertexEnumerator();
                    output = enumerator.enumerate(working);
                    stats = enumerator.getLastStats();
                    // Write output in lrslib-compatible format
                    java.io.PrintWriter pw = new java.io.PrintWriter(System.out, true);
                    output.write(pw);
                    pw.flush();
                    break;
                }
                case CH: { // V -> H (facet enumeration)
                    FacetEnumerator enumerator = new FacetEnumerator();
                    output = enumerator.enumerate(working);
                    stats = enumerator.getLastStats();
                    try (java.io.PrintWriter pw = new java.io.PrintWriter(System.out)) {
                        output.write(pw);
                    }
                    break;
                }
                default:
                    throw new IllegalStateException("Unknown mode: " + dat.mode);
            }

            final Instant t1 = Instant.now();
            double secs = Duration.between(t0, t1).toMillis() / 1000.0;

            // 4) Print stats (mirrors lrslib totals; adjust fields as you add parity)
            // After writing the output polyhedron
            if (stats != null) {
                if (dat.mode == LrsDat.Mode.VE) {
                    stats.printVertexTotals();
                }
                if (dat.printCobasis && stats.lastCobasis != null) {
                    System.out.print("printcobasis 1\n");
                    System.out.println(java.util.Arrays.toString(stats.lastCobasis));
                }
            }
            System.out.printf("*elapsed time: %.3f seconds%n", secs);

            return 0;
        } catch (FileNotFoundException e) {
            System.err.println("File not found: " + inputPath);
            return 1;
        } catch (IOException e) {
            System.err.println("I/O error: " + e.getMessage());
            return 1;
        } catch (RuntimeException e) {
            // Use this block to map certain exceptions to a restart in the future
            // (e.g., resource limits, user-requested checkpointing). For now, just report.
            System.err.println("*unrecoverable error: " + e.getMessage());
            return -1;
        }
    }

    /**
     * Utility that provides the same help text as your current Main. Kept separate
     * so LrsDriver can be used standalone or from Main.
     */
    static final class MainErrorPrinter {
        private static void usage() {
            System.err.println(
                    "Usage: lrs-java [options] <input-file>\n" +
                            "Modes:\n" +
                            "  -v            H->V (vertex enumeration)  [default]\n" +
                            "  -h            V->H (facet enumeration)\n" +
                            "Options:\n" +
                            "  -redund       remove redundant inequalities/generators (if implemented)\n" +
                            "  -minrep       minimize representation (hidden linearities, etc.)\n" +
                            "  -printcobasis print the final cobasis\n" +
                            "  -integer      declare input as integer data (affects stats only)\n" +
                            "  -seed <n>     random seed (if you later add randomized tie-breaks)\n" +
                            "  -maxdepth <d> limit reverse-search depth (debug)\n" +
                            "  -threads <t>  parallel facet/vertex post-processing (future)\n" +
                            "  -eliminate i,j,k  eliminate the listed columns (projection)\n" +
                            "  -linset i,j,k    mark listed constraint rows as equalities\n" +
                            "  -project i,j,k   project to these columns (kept set)\n"
            );
        }
    }
}

/**
 * Java analogue of lrs_restart_dat (lrsrestart.h) and the bookkeeping that lives in lrsdriver.c.
 * Only fields that make sense for Java are included; feel free to extend as you bring over
 * more of lrslib's counters and options.
 */
final class RestartData {
    // Behavior flags
    boolean overrideQ = false; // overide in C (spelled that way in lrslib)
    boolean restart = false;   // whether to restart a run (unused in BigInteger arithmetic)

    // Dimensions & limits
    int d = 0;                 // current dimension (columns excluding homogenizing constant)
    int maxCobases = 0;        // optional cap on cobases visited
    int maxDepth = -1;         // search depth limit (-1 means unlimited)
    int minDepth = 0;          // minimum depth to start collecting items

    // Counters (10 slots mirror lrsdriver usage; repurpose as needed)
    final long[] count = new long[10];

    // Live traversal state
    int depth = 0;

    // Modes/verbosity
    int lrs = 1;               // 1 means enable lrs core
    int redund = 0;            // redundancy pass on/off
    int messages = 0;          // verbosity level

    // Additional switches carried around by lrslib
    int fel = 0;               // facet enumeration level (if used)
    int testlin = 0;           // test linearity
    int redundphase = 0;       // which redundancy phase

    // Redundancy helpers (row marks, etc.)
    boolean[] redineq = null;

    void resetDefaults() {
        overrideQ = false;
        restart = false;
        d = 0;
        maxCobases = 0;
        maxDepth = -1;
        minDepth = 0;
        Arrays.fill(count, 0L);
        depth = 0;
        lrs = 1;
        redund = 0;
        messages = 0;
        fel = 0;
        testlin = 0;
        redundphase = 0;
        redineq = null;
    }
}
