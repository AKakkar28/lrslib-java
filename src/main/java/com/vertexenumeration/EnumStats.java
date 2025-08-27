package com.vertexenumeration;

public final class EnumStats {
    public int vertices;
    public int rays;
    public int bases;
    public int integerVertices;
    public int maxDepth;
    public int[] lastCobasis;

    public void printVertexTotals() {
        System.out.printf("*Totals: vertices=%d rays=%d bases=%d integer_vertices=%d%n",
                vertices, rays, bases, integerVertices);
        System.out.printf("*max_vertex_depth=%d%n", maxDepth);
    }
}
