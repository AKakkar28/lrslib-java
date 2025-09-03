package com.vertexenumeration;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;
import java.util.*;

public class SimplexDictionaryTest {

    /** Unit square in R^2:
     *  H rows are b + A x >= 0 for:
     *   x >= 0, y >= 0, 1 - x >= 0, 1 - y >= 0
     *   → rows: [0,1,0], [0,0,1], [1,-1,0], [1,0,-1]
     */
    private static Fraction[][] unitSquareH() {
        Fraction Z = Fraction.ZERO, O = Fraction.ONE;
        return new Fraction[][]{
                { Z,  O,  Z },
                { Z,  Z,  O },
                { O,  Fraction.of(-1),  Z },
                { O,  Z,  Fraction.of(-1) }
        };
    }



}
