package com.vertexenumeration;

import java.util.*;

public final class LrsDat {
    public enum Mode { VE, CH }             // vertices / facets (convex hull)
    public final Mode mode;

    // flags / features
    public final boolean redund;            // remove redundant inequalities
    public final boolean minrep;            // hidden linearities + minimal H
    public final boolean printCobasis;      // print cobasis alongside output
    public final int maxDepth;              // reverse-search depth cap (0 = no cap)
    public final int threads;               // root parallelism only (like 7.x)
    public final long seed;                 // tie-breaking, deterministic

    // transforms
    public final int[] eliminate;           // column indices to eliminate (fel)
    public final int[] project;             // columns to keep (project)
    public final int[] linearities;         // indices of linear (==) rows
    public final boolean integerInput;      // input declared integer

    private LrsDat(Builder b) {
        this.mode = b.mode;
        this.redund = b.redund;
        this.minrep = b.minrep;
        this.printCobasis = b.printCobasis;
        this.maxDepth = b.maxDepth;
        this.threads = b.threads;
        this.seed = b.seed;
        this.eliminate = b.eliminate.toArray();
        this.project = b.project.toArray();
        this.linearities = b.linearities.toArray();
        this.integerInput = b.integerInput;
    }

    public static final class Builder {
        private Mode mode = Mode.VE;
        private boolean redund, minrep, printCobasis, integerInput;
        private int maxDepth = 0, threads = 1;
        private long seed = 1L;
        private final IntSet eliminate = new IntSet();
        private final IntSet project = new IntSet();
        private final IntSet linearities = new IntSet();

        public Builder mode(Mode m){ this.mode=m; return this; }
        public Builder redund(boolean v){ this.redund=v; return this; }
        public Builder minrep(boolean v){ this.minrep=v; return this; }
        public Builder printCobasis(boolean v){ this.printCobasis=v; return this; }
        public Builder maxDepth(int v){ this.maxDepth=v; return this; }
        public Builder threads(int v){ this.threads=Math.max(1,v); return this; }
        public Builder seed(long v){ this.seed=v; return this; }
        public Builder integerInput(boolean v){ this.integerInput=v; return this; }
        public Builder addEliminate(int c){ this.eliminate.add(c); return this; }
        public Builder addProject(int c){ this.project.add(c); return this; }
        public Builder addLinearity(int r){ this.linearities.add(r); return this; }
        public LrsDat build(){ return new LrsDat(this); }
    }

    // tiny sorted, unique int bag
    private static final class IntSet {
        private final TreeSet<Integer> s = new TreeSet<>();
        void add(int v){ s.add(v); }
        int[] toArray(){ int[] a=new int[s.size()]; int i=0; for(int v:s)a[i++]=v; return a; }
    }
}
