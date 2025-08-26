package com.vertexenumeration;

import static org.junit.jupiter.api.Assertions.*;

import java.math.BigInteger;

import org.junit.jupiter.api.Test;

/**
 * Unit tests for the {@link Fraction} class.  These tests verify basic
 * arithmetic operations, equality, and comparison.
 */
public class FractionTest {

    @Test
    public void testArithmetic() {
        Fraction a = new Fraction(new BigInteger("1"), new BigInteger("2"));
        Fraction b = new Fraction(new BigInteger("3"), new BigInteger("4"));
        // 1/2 + 3/4 = 5/4
        assertEquals(new Fraction(new BigInteger("5"), new BigInteger("4")), a.add(b));
        // 1/2 - 3/4 = -1/4
        assertEquals(new Fraction(new BigInteger("-1"), new BigInteger("4")), a.subtract(b));
        // 1/2 * 3/4 = 3/8
        assertEquals(new Fraction(new BigInteger("3"), new BigInteger("8")), a.multiply(b));
        // (1/2) / (3/4) = 2/3
        assertEquals(new Fraction(new BigInteger("2"), new BigInteger("3")), a.divide(b));
    }

    @Test
    public void testEqualityAndComparison() {
        Fraction x = new Fraction(new BigInteger("2"), new BigInteger("4"));
        Fraction y = new Fraction(new BigInteger("1"), new BigInteger("2"));
        assertEquals(x, y);
        assertEquals(0, x.compareTo(y));
        Fraction z = new Fraction(new BigInteger("3"), new BigInteger("2"));
        assertTrue(z.compareTo(y) > 0);
        assertTrue(y.compareTo(z) < 0);
    }
}