package com.vertexenumeration;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

/** Tests for the tableau-based Dictionary (pivot + lex ratio rule). */
public class DictionaryTest {

    /** Build a tiny feasible dictionary with two constraints & two nonbasics.
     *  Row/col indexing matches Dictionary docs: row0=objective, col0=RHS. */
    private static Dictionary tinyLexTableau() {
        // Reduced costs in row0 are negative → we must improve
        // T = [
        //   [ 0 | -1  -1 ]
        //   [ 1 |  1   0 ]
        //   [ 1 |  0   1 ]
        // ]
        Fraction Z = Fraction.ZERO, O = Fraction.ONE;
        Fraction[][] T = {
                { Z, Fraction.of(-1), Fraction.of(-1) },
                { O, O, Z },
                { O, Z, O }
        };
        // basis rows (size m+1, basis[0] unused): slacks initially basic
        int[] basis   = { 0, 3, 4 };
        // cobasis cols (size n+1, cobasis[0] unused): original variables 1..n
        int[] cobasis = { 0, 1, 2 };
        Dictionary D = Dictionary.of(T, basis, cobasis);
        D.setPivotRule(Dictionary.PivotRule.LEXICOGRAPHIC);
        return D;
    }

    @Test
    public void testFeasibilityAndNegativeReducedCostDetection() {
        Dictionary D = tinyLexTableau();
        assertTrue(D.isFeasible(), "RHS >= 0 ⇒ feasible");
        assertEquals(Dictionary.LPStatus.RUNNING, D.status());
    }

    @Test
    public void testSinglePivotAndLabelSwap() {
        Dictionary D = tinyLexTableau();

        // Entering: first col with negative reduced cost → col 1
        int enter = 1;
        // Leaving (lex rule): only row1 has positive a_re (1>0) → row1
        int leave = D.leavingRowFor(enter);
        assertEquals(1, leave);

        // Perform pivot and check labels swap correctly
        D.pivot(leave, enter);
        int[] basis = D.basis();
        int[] cob   = D.cobasis();
        assertEquals(1, basis[1], "var x1 enters basis");
        assertEquals(3, cob[1],   "slack s1 leaves basis");

        // Tableau should now have pivot col as unit vector and pivot row RHS normalized
        Fraction[][] T = D.tableau();
        assertEquals(Fraction.ONE,  T[1][1]);
        assertEquals(Fraction.ZERO, T[0][1]);
        assertEquals(Fraction.ZERO, T[2][1]);
    }

    @Test
    public void testSolveToOptimal() {
        Dictionary D = tinyLexTableau();
        Dictionary.LPStatus st = D.solve(20);
        assertEquals(Dictionary.LPStatus.OPTIMAL, st);
    }

    @Test
    public void testUnboundedDetection() {
        // Make both reduced costs negative but column 2 nonpositive in all rows → unbounded
        Fraction Z = Fraction.ZERO, O = Fraction.ONE;
        Fraction[][] T = {
                { Z, Fraction.of(-1), Fraction.of(-1) },
                { O, O, Z },   // col2 = 0
                { O, Z, Z }    // col2 = 0 ⇒ no leaving row for enter=2
        };
        int[] basis   = { 0, 3, 4 };
        int[] cobasis = { 0, 1, 2 };
        Dictionary D = Dictionary.of(T, basis, cobasis);
        D.setPivotRule(Dictionary.PivotRule.LEXICOGRAPHIC);
        D.step(); // should pick enter=1 (since -1 at col1), not unbounded yet
        // Force enter=2 by zeroing row0 col1
        D.tableau()[0][1] = Z; // NOTE: if tableau() returns a copy in your impl, instead rebuild another D
        // Rebuild with col1 cost zeroed to drive enter=2 scenario:
        Fraction[][] T2 = {
                { Z, Z, Fraction.of(-1) },
                { O, O, Z },
                { O, Z, Z }
        };
        D = Dictionary.of(T2, basis, cobasis);
        D.setPivotRule(Dictionary.PivotRule.LEXICOGRAPHIC);
        D.step();
        assertEquals(Dictionary.LPStatus.UNBOUNDED, D.status());
    }

    @Test
    public void testBlandVsLexTieBreak() {
        // Degenerate tie on ratios: both rows give same ratio; Bland chooses smaller row index.
        Fraction Z = Fraction.ZERO, O = Fraction.ONE;
        Fraction[][] T = {
                { Z, Fraction.of(-1), Z },    // enter col1
                { O, Fraction.of(1), Z },
                { O, Fraction.of(1), Z }      // identical row → tie
        };
        int[] basis   = { 0, 3, 4 };
        int[] cobasis = { 0, 1, 2 };

        Dictionary bland = Dictionary.of(T, basis, cobasis);
        bland.setPivotRule(Dictionary.PivotRule.BLAND);
        assertEquals(1, bland.leavingRowFor(1), "Bland → smallest row index");

        Dictionary lex = Dictionary.of(T, basis, cobasis);
        lex.setPivotRule(Dictionary.PivotRule.LEXICOGRAPHIC);
        assertEquals(1, lex.leavingRowFor(1), "Lex also falls back to smaller row on perfect tie");
    }
}
