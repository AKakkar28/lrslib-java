package com.vertexenumeration;

import java.util.Arrays;

/** Snapshot of reverse-search state at a node (basis + last move). */
final class RSState {
    final int[] basis;     // basis labels for rows 1..m (basis[0] ignored)
    final int[] cobasis;   // cobasis labels for cols 1..n (cobasis[0] ignored)
    final int lastLeaveRow;
    final int lastEnterCol;

    RSState(int[] basis, int[] cobasis, int lastLeaveRow, int lastEnterCol) {
        this.basis = Arrays.copyOf(basis, basis.length);
        this.cobasis = Arrays.copyOf(cobasis, cobasis.length);
        this.lastLeaveRow = lastLeaveRow;
        this.lastEnterCol = lastEnterCol;
    }
}
