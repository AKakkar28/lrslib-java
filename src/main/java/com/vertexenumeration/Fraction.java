package com.vertexenumeration;

import java.math.BigInteger;
import java.util.Objects;

/**
 * Immutable arbitrary‑precision rational number.  Internally stores a
 * numerator and denominator as {@link BigInteger}s and reduces to
 * lowest terms on construction.  Denominator is always positive.
 */
public final class Fraction implements Comparable<Fraction> {
    private final BigInteger numerator;
    private final BigInteger denominator;

    /**
     * Constructs a {@code Fraction} given numerator and denominator.  The
     * fraction is reduced and the denominator is forced to be positive.
     *
     * @param num the numerator
     * @param den the denominator (must not be zero)
     */
    public Fraction(BigInteger num, BigInteger den) {
        Objects.requireNonNull(num, "numerator");
        Objects.requireNonNull(den, "denominator");
        if (den.equals(BigInteger.ZERO)) {
            throw new ArithmeticException("Zero denominator");
        }
        // Normalize sign: denominator always positive
        if (den.signum() < 0) {
            num = num.negate();
            den = den.negate();
        }
        // Reduce fraction
        BigInteger gcd = num.gcd(den);
        this.numerator = num.divide(gcd);
        this.denominator = den.divide(gcd);
    }

    /** Constructs a {@code Fraction} representing an integer. */
    public Fraction(BigInteger integer) {
        this(integer, BigInteger.ONE);
    }

    /** Returns the numerator. */
    public BigInteger getNumerator() {
        return numerator;
    }

    /** Returns the denominator. */
    public BigInteger getDenominator() {
        return denominator;
    }

    /** Returns a new {@code Fraction} equal to this plus {@code other}. */
    public Fraction add(Fraction other) {
        BigInteger num = this.numerator.multiply(other.denominator)
                                .add(other.numerator.multiply(this.denominator));
        BigInteger den = this.denominator.multiply(other.denominator);
        return new Fraction(num, den);
    }

    /** Returns a new {@code Fraction} equal to this minus {@code other}. */
    public Fraction subtract(Fraction other) {
        BigInteger num = this.numerator.multiply(other.denominator)
                                .subtract(other.numerator.multiply(this.denominator));
        BigInteger den = this.denominator.multiply(other.denominator);
        return new Fraction(num, den);
    }

    /** Returns a new {@code Fraction} equal to this multiplied by {@code other}. */
    public Fraction multiply(Fraction other) {
        BigInteger num = this.numerator.multiply(other.numerator);
        BigInteger den = this.denominator.multiply(other.denominator);
        return new Fraction(num, den);
    }

    /** Returns a new {@code Fraction} equal to this divided by {@code other}. */
    public Fraction divide(Fraction other) {
        if (other.numerator.equals(BigInteger.ZERO)) {
            throw new ArithmeticException("Divide by zero fraction");
        }
        BigInteger num = this.numerator.multiply(other.denominator);
        BigInteger den = this.denominator.multiply(other.numerator);
        return new Fraction(num, den);
    }

    @Override
    public int compareTo(Fraction other) {
        // Compare by cross multiplication to avoid rounding
        BigInteger lhs = this.numerator.multiply(other.denominator);
        BigInteger rhs = other.numerator.multiply(this.denominator);
        return lhs.compareTo(rhs);
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof Fraction)) return false;
        Fraction o = (Fraction) obj;
        return numerator.equals(o.numerator) && denominator.equals(o.denominator);
    }

    @Override
    public int hashCode() {
        return Objects.hash(numerator, denominator);
    }

    @Override
    public String toString() {
        if (denominator.equals(BigInteger.ONE)) {
            return numerator.toString();
        }
        return numerator + "/" + denominator;
    }
}