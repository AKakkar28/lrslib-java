package com.vertexenumeration;

import java.math.BigInteger;
import java.util.Objects;

/** Immutable arbitrary-precision rational with normalized sign and gcd reduction. */
public final class Fraction implements Comparable<Fraction> {
    public static final Fraction ZERO = new Fraction(BigInteger.ZERO, BigInteger.ONE);
    public static final Fraction ONE  = new Fraction(BigInteger.ONE,  BigInteger.ONE);

    private final BigInteger numerator;
    private final BigInteger denominator; // always > 0

    public Fraction(BigInteger num, BigInteger den) {
        Objects.requireNonNull(num, "numerator");
        Objects.requireNonNull(den, "denominator");
        if (den.signum() == 0) throw new ArithmeticException("Zero denominator");
        if (den.signum() < 0) { num = num.negate(); den = den.negate(); }
        BigInteger g = num.gcd(den);
        this.numerator = num.divide(g);
        this.denominator = den.divide(g);
    }

    public Fraction(BigInteger integer) { this(integer, BigInteger.ONE); }

    /** Convenience factories */
    public static Fraction of(long n) { return new Fraction(BigInteger.valueOf(n), BigInteger.ONE); }
    public static Fraction of(long n, long d) { return new Fraction(BigInteger.valueOf(n), BigInteger.valueOf(d)); }

    public BigInteger getNumerator()   { return numerator; }
    public BigInteger getDenominator() { return denominator; }

    public Fraction add(Fraction o) {
        BigInteger num = numerator.multiply(o.denominator).add(o.numerator.multiply(denominator));
        BigInteger den = denominator.multiply(o.denominator);
        return new Fraction(num, den);
    }
    public Fraction subtract(Fraction o) {
        BigInteger num = numerator.multiply(o.denominator).subtract(o.numerator.multiply(denominator));
        BigInteger den = denominator.multiply(o.denominator);
        return new Fraction(num, den);
    }
    public Fraction multiply(Fraction o) {
        return new Fraction(numerator.multiply(o.numerator), denominator.multiply(o.denominator));
    }
    public Fraction divide(Fraction o) {
        if (o.numerator.signum() == 0) throw new ArithmeticException("Divide by zero fraction");
        return new Fraction(numerator.multiply(o.denominator), denominator.multiply(o.numerator));
    }

    public Fraction negate() { return numerator.signum() == 0 ? ZERO : new Fraction(numerator.negate(), denominator); }
    public Fraction abs()    { return numerator.signum() < 0 ? negate() : this; }
    public int signum()      { return numerator.signum(); }
    public boolean isZero()  { return numerator.signum() == 0; }

    @Override public int compareTo(Fraction o) {
        // cross-multiply: a/b ? c/d  <=>  ad ? cb
        return numerator.multiply(o.denominator).compareTo(o.numerator.multiply(denominator));
    }
    @Override public boolean equals(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof Fraction)) return false;
        Fraction o = (Fraction) obj;
        return numerator.equals(o.numerator) && denominator.equals(o.denominator);
    }
    @Override public int hashCode() { return Objects.hash(numerator, denominator); }
    @Override public String toString() {
        return denominator.equals(BigInteger.ONE) ? numerator.toString() : numerator + "/" + denominator;
    }
}
