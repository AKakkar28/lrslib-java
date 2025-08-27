package com.vertexenumeration;

public interface LrsDic {
    int m();                   // constraints
    int n();                   // nonbasic vars
    Fraction[][] tableau();    // copy or view
    int[] basis();             // labels of basic variables
    int[] cobasis();           // labels of nonbasic variables
    boolean isFeasible();
    void setPivotRule(Dictionary.PivotRule rule);
    int leavingRowFor(int enterCol);
    boolean canPivot(int enterCol);
    void pivot(int leaveRow, int enterCol);
    Dictionary.LPStatus solve(int maxIters);
    Dictionary.LPStatus status();
}
