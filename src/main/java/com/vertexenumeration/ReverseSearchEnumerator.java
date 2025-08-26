package com.vertexenumeration;

import java.util.*;

final class ReverseSearchEnumerator {

    static final class Result {
        final List<Fraction[]> vertices;
        final int bases;
        final int maxDepth;
        Result(List<Fraction[]> v, int b, int d) { vertices = v; bases = b; maxDepth = d; }
    }

    static final class Basis {
        final int[] rows;      // sorted tight-constraint indices (size d)
        final Fraction[] x;    // vertex coords (length d), solves A x = -a0
        Basis(int[] rows, Fraction[] x) { this.rows = rows; this.x = x; }
        int compareTo(Basis o) { // lex compare rows
            for (int i = 0; i < rows.length; i++) {
                int c = Integer.compare(this.rows[i], o.rows[i]);
                if (c != 0) return c;
            }
            return 0;
        }
    }

    private final Fraction[][] H; // m x n (n = 1 + d)
    private final int m, n, d;
    private final Fraction ZERO, ONE;

    ReverseSearchEnumerator(Fraction[][] H) {
        this.H = H;
        this.m = H.length;
        this.n = (m == 0 ? 0 : H[0].length);
        this.d = Math.max(0, n - 1);
        if (m == 0 || n == 0 || d == 0) throw new IllegalArgumentException("Empty/degenerate matrix");
        this.ZERO = anyZero(H);
        this.ONE  = anyOne(H, ZERO);
    }

    /** Enumerate vertices via reverse search over feasible bases; return vertices + stats. */
    Result run() {
        // 1) gather all feasible bases (nondegenerate vertices)
        List<Basis> bases = new ArrayList<>();
        int[] comb = initComb(d);
        while (comb != null) {
            Basis B = solveIfFeasible(comb);
            if (B != null) bases.add(B);
            comb = nextComb(comb, m, d);
        }
        if (bases.isEmpty()) return new Result(Collections.emptyList(), 0, 0);

        Map<String,Basis> basisMap = new HashMap<>(bases.size()*2);
        for (Basis b : bases) basisMap.put(key(b.rows), b);

        // 2) choose root as lex-min feasible basis
        bases.sort(Basis::compareTo);
        Basis root = bases.get(0);

        // parent pointers computed lazily; root has none
        Map<String,String> parent = new HashMap<>();
        parent.put(key(root.rows), "ROOT");

        List<Fraction[]> out = new ArrayList<>();
        int basesVisited = 0, maxDepth = 0;

        // 3) DFS with explicit depth, visiting children in strict lex order
        Deque<Basis> stack = new ArrayDeque<>();
        Deque<Integer> depth = new ArrayDeque<>();
        stack.push(root); depth.push(0);

        while (!stack.isEmpty()) {
            Basis cur = stack.pop();
            int dep = depth.pop();
            basesVisited++;
            if (dep > maxDepth) maxDepth = dep;

            // output vertex (homogeneous [1, x...])
            out.add(toHomogeneous(cur.x));

            // neighbors whose parent is this basis, visited in lex order
            TreeSet<Basis> children = new TreeSet<>(Basis::compareTo);
            BitSet inBasis = bitsetOf(cur.rows, m);

            for (int r = 0; r < m; r++) if (!inBasis.get(r)) {
                for (int pos = 0; pos < d; pos++) {
                    int[] neighRows = cur.rows.clone();
                    neighRows[pos] = r;
                    Arrays.sort(neighRows);
                    Basis neigh = basisMap.get(key(neighRows));
                    if (neigh == null || neigh.compareTo(cur) == 0) continue;
                    String k = key(neigh.rows);
                    if (!parent.containsKey(k)) {
                        Basis p = computeParent(neigh, basisMap);
                        if (p != null && key(p.rows).equals(key(cur.rows))) {
                            parent.put(k, key(cur.rows));
                            children.add(neigh);
                        }
                    }
                }
            }

            // push in reverse so smallest lex child is popped first
            for (Basis child : children.descendingSet()) {
                stack.push(child); depth.push(dep + 1);
            }
        }

        return new Result(out, basesVisited, maxDepth);
    }

    // -------- internals --------

    private Basis solveIfFeasible(int[] rows) {
        Fraction[][] A = new Fraction[d][d];
        Fraction[] b = new Fraction[d];
        for (int i = 0; i < d; i++) {
            int r = rows[i];
            b[i] = ZERO.subtract(H[r][0]);  // -a0
            for (int j = 0; j < d; j++) A[i][j] = H[r][j+1];
        }
        Fraction[] x = solve(A, b, ZERO);
        if (x == null || !feasible(H, x, ZERO)) return null;
        return new Basis(sortedCopy(rows), x);
    }

    /** Parent(neigh) = lexicographically smallest neighbor that is strictly < neigh. */
    private Basis computeParent(Basis neigh, Map<String,Basis> basisMap) {
        BitSet in = bitsetOf(neigh.rows, m);
        Basis best = null;
        for (int r = 0; r < m; r++) if (!in.get(r)) {
            for (int pos = 0; pos < d; pos++) {
                int[] candRows = neigh.rows.clone();
                candRows[pos] = r;
                Arrays.sort(candRows);
                Basis cand = basisMap.get(key(candRows));
                if (cand == null) continue;
                if (cand.compareTo(neigh) < 0 && (best == null || cand.compareTo(best) < 0)) {
                    best = cand;
                }
            }
        }
        return best;
    }

    private Fraction[] toHomogeneous(Fraction[] x) {
        Fraction[] v = new Fraction[d+1];
        v[0] = ONE;
        for (int j = 0; j < d; j++) v[j+1] = x[j];
        return v;
    }

    private static BitSet bitsetOf(int[] rows, int m) { BitSet bs = new BitSet(m); for (int r: rows) bs.set(r); return bs; }
    private static int[] sortedCopy(int[] a) { int[] b = a.clone(); Arrays.sort(b); return b; }
    private static String key(int[] rows) { StringBuilder sb=new StringBuilder(rows.length*3); for(int r:rows) sb.append(r).append(','); return sb.toString(); }

    private static Fraction anyZero(Fraction[][] H) { for (Fraction[] row : H) for (Fraction f : row) return f.subtract(f); throw new IllegalStateException("Empty H"); }
    private static Fraction anyOne (Fraction[][] H, Fraction ZERO) { for (Fraction[] row : H) for (Fraction f : row) if (f.compareTo(ZERO)!=0) return f.divide(f); throw new IllegalStateException("Cannot build ONE"); }

    private static boolean feasible(Fraction[][] H, Fraction[] x, Fraction ZERO) {
        for (Fraction[] row : H) {
            Fraction lhs = row[0];
            for (int j = 0; j < x.length; j++) lhs = lhs.add(row[j+1].multiply(x[j]));
            if (lhs.compareTo(ZERO) < 0) return false;
        }
        return true;
    }

    private static Fraction[] solve(Fraction[][] A, Fraction[] b, Fraction ZERO) {
        int n = b.length;
        Fraction[][] M = new Fraction[n][n+1];
        for (int i = 0; i < n; i++) { System.arraycopy(A[i], 0, M[i], 0, n); M[i][n] = b[i]; }
        int row = 0;
        for (int col = 0; col < n && row < n; col++) {
            int piv = row; while (piv < n && M[piv][col].compareTo(ZERO) == 0) piv++;
            if (piv == n) continue;
            if (piv != row) { Fraction[] t = M[piv]; M[piv]=M[row]; M[row]=t; }
            Fraction diag = M[row][col]; if (diag.compareTo(ZERO) == 0) continue;
            for (int j = col; j <= n; j++) M[row][j] = M[row][j].divide(diag);
            for (int r = 0; r < n; r++) if (r != row) {
                Fraction factor = M[r][col];
                if (factor.compareTo(ZERO) != 0) for (int j = col; j <= n; j++) M[r][j] = M[r][j].subtract(factor.multiply(M[row][j]));
            }
            row++;
        }
        for (int i = 0; i < n; i++) { boolean allZero = true; for (int j=0;j<n;j++) if (M[i][j].compareTo(ZERO)!=0){allZero=false;break;} if (allZero && M[i][n].compareTo(ZERO)!=0) return null; }
        Fraction[] x = new Fraction[n];
        for (int i = 0; i < n; i++) { int lead = -1; for (int j=0;j<n;j++) if (M[i][j].compareTo(ZERO)!=0){lead=j;break;} if (lead==-1) return null; x[lead] = M[i][n]; }
        for (int i = 0; i < n; i++) if (x[i] == null) return null;
        return x;
    }

    private static int[] initComb(int k){ int[] a=new int[k]; for(int i=0;i<k;i++) a[i]=i; return a; }
    private static int[] nextComb(int[] a,int n,int k){ int i=k-1; while(i>=0 && a[i]==n-k+i) i--; if(i<0) return null; a[i]++; for(int j=i+1;j<k;j++) a[j]=a[j-1]+1; return a; }
}
