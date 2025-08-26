package com.vertexenumeration;

public final class EnumStats {
    public int vertices;
    public int rays;
    public int bases;
    public int integerVertices;
    public int maxDepth;

    @Override public String toString() {
        return "*Totals: vertices=" + vertices +
                " rays=" + rays +
                " bases=" + bases +
                " integer_vertices=" + integerVertices +
                "  max_vertex_depth=" + maxDepth;
    }
}
