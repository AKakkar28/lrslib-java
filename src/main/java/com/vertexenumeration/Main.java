package com.vertexenumeration;

public final class Main {

    // Keep the usage text here for convenience; LrsDriver handles parsing/execution.
    @SuppressWarnings("unused")
    private static void usage() {
        System.err.println(
                "Usage: lrs-java [options] <input-file>\n" +
                        "Modes:\n" +
                        "  -v            H->V (vertex enumeration)  [default]\n" +
                        "  -h            V->H (facet enumeration)\n" +
                        "Options:\n" +
                        "  -redund                 remove redundant inequalities (H only)\n" +
                        "  -minrep                 hidden linearities + minimal H (H only)\n" +
                        "  -eliminate c1,c2,...    Fourier–Motzkin eliminate columns (H only)\n" +
                        "  -project   c1,c2,...    keep only these columns (H only)\n" +
                        "  -linearity r1,r2,...    mark rows as equalities (H only)\n" +
                        "  -printcobasis           print cobasis indices with output\n" +
                        "  -seed N                 deterministic tie-breaking seed\n" +
                        "  -maxdepth N             cap reverse-search depth (0 = no cap)\n" +
                        "  -threads N              root parallelism (implementation-defined)\n" +
                        "  -integer                declare input as integer data (affects stats only)\n"
        );
    }

    public static void main(String[] args) {
        int code = new LrsDriver().run(args);
        System.exit(code);
    }
}
