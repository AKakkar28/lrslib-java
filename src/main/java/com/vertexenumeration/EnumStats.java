package com.vertexenumeration;

public final class EnumStats {
    public int vertices;
    public int rays;
    public int bases;
    public int facets;
    public int integerVertices;
    public int maxDepth;
    public int[] lastCobasis;

    private LrsDat.Mode mode;

    public void setMode(LrsDat.Mode m) {
        this.mode = m;
    }

    @Override
    public String toString() {
        if (mode == LrsDat.Mode.VE) {
            return "*Totals: vertices=" + vertices +
                    " rays=" + rays +
                    " bases=" + bases +
                    " integer_vertices=" + integerVertices +
                    "  max_vertex_depth=" + maxDepth;
        } else if (mode == LrsDat.Mode.CH) {
            return "*Totals: facets=" + facets +
                    " bases=" + bases;
        } else {
            return "*Totals: (unknown mode)";
        }
    }
}
