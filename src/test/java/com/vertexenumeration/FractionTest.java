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

    @Test
    public void testNormalizationAndSign() {
        // denominator always positive; gcd reduced
        assertEquals(new Fraction(BigInteger.valueOf(-2), BigInteger.valueOf(3)),
                new Fraction(BigInteger.valueOf(4),  BigInteger.valueOf(-6)));
        assertEquals(new Fraction(BigInteger.ONE, BigInteger.ONE),
                new Fraction(BigInteger.valueOf(10), BigInteger.valueOf(10)));
        // zero normalizes to 0/1
        assertEquals(Fraction.ZERO, new Fraction(BigInteger.ZERO, BigInteger.valueOf(-7)));
        assertEquals(BigInteger.ONE, Fraction.ZERO.denominator());
    }

    @Test
    public void testZeroOneConstants() {
        assertEquals(0, Fraction.ZERO.signum());
        assertTrue(Fraction.ZERO.isZero());
        assertEquals(1, Fraction.ONE.signum());
        assertFalse(Fraction.ONE.isZero());
        assertEquals("0", Fraction.ZERO.toString());
        assertEquals("1", Fraction.ONE.toString());
    }

    @Test
    public void testNegateAbsSignum() {
        Fraction p = new Fraction(BigInteger.valueOf(5), BigInteger.valueOf(7));
        Fraction n = p.negate();
        assertEquals(-1, n.signum());
        assertEquals(p, n.abs());
        assertEquals(1, p.signum());
        assertEquals(0, Fraction.ZERO.negate().signum());
    }

    @Test
    public void testParseAndToString() {
        assertEquals(new Fraction(BigInteger.valueOf(3), BigInteger.valueOf(4)),
                Fraction.parse(" 3/4 "));
        assertEquals(new Fraction(BigInteger.valueOf(-2), BigInteger.ONE),
                Fraction.parse("-2"));
        assertEquals("7/3", new Fraction(BigInteger.valueOf(7), BigInteger.valueOf(3)).toString());
        assertEquals("5",   new Fraction(BigInteger.valueOf(5), BigInteger.ONE).toString());
    }

    @Test
    public void testDivideByZeroThrows() {
        Fraction a = new Fraction(BigInteger.ONE, BigInteger.TWO);
        assertThrows(ArithmeticException.class, () -> a.divide(Fraction.ZERO));
    }

    @Test
    public void testInverse() {
        Fraction f = new Fraction(BigInteger.valueOf(-3), BigInteger.valueOf(7));
        assertEquals(new Fraction(BigInteger.valueOf(-7), BigInteger.valueOf(3)), f.inverse());
        assertThrows(ArithmeticException.class, () -> Fraction.ZERO.inverse());
    }

    @Test
    public void testImmutability() {
        Fraction a = new Fraction(BigInteger.ONE, BigInteger.valueOf(3));
        Fraction b = new Fraction(BigInteger.ONE, BigInteger.valueOf(6));
        Fraction sum = a.add(b);          // == 1/2
        assertEquals(new Fraction(BigInteger.ONE, BigInteger.valueOf(3)), a); // a unchanged
        assertEquals(new Fraction(BigInteger.ONE, BigInteger.valueOf(6)), b); // b unchanged
        assertEquals(new Fraction(BigInteger.ONE, BigInteger.valueOf(2)), sum);
    }

    @Test
    public void testComparisonTransitivityAndTotalOrder() {
        Fraction a = Fraction.parse("1/3");
        Fraction b = Fraction.parse("2/5");
        Fraction c = Fraction.parse("3/7");
        // a < b < c ?
        assertTrue(a.compareTo(b) < 0);
        assertTrue(b.compareTo(c) < 0);
        assertTrue(a.compareTo(c) < 0); // transitive
        // equals implies compareTo == 0
        assertEquals(0, a.compareTo(Fraction.parse("2/6")));
    }


    @Test
    public void testHashEqualsContract() {
        Fraction a = new Fraction(BigInteger.valueOf(10), BigInteger.valueOf(20));
        Fraction b = new Fraction(BigInteger.ONE, BigInteger.valueOf(2));
        assertEquals(a, b);
        assertEquals(a.hashCode(), b.hashCode());
    }



}