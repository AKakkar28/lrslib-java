package com.vertexenumeration;

import java.math.BigInteger;
import java.util.Objects;

/** Immutable arbitrary-precision rational with normalized sign and gcd reduction. */
public final class Fraction implements Numeric<Fraction> {
    public static final Fraction ZERO = new Fraction(BigInteger.ZERO, BigInteger.ONE);
    public static final Fraction ONE  = new Fraction(BigInteger.ONE,  BigInteger.ONE);

    private final BigInteger n;        // numerator
    private final BigInteger d;        // denominator > 0

    /** Creates and reduces; denominator must be nonzero. */
    public Fraction(BigInteger num, BigInteger den) {
        Objects.requireNonNull(num, "numerator");
        Objects.requireNonNull(den, "denominator");
        if (den.signum() == 0) throw new ArithmeticException("Zero denominator");
        if (den.signum() < 0) { num = num.negate(); den = den.negate(); }
        BigInteger g = num.gcd(den);
        this.n = num.divide(g);
        this.d = den.divide(g);
    }

    public Fraction(BigInteger integer) { this(integer, BigInteger.ONE); }

    /** Factories */
    public static Fraction of(long k) { return new Fraction(BigInteger.valueOf(k), BigInteger.ONE); }
    public static Fraction of(long num, long den) { return new Fraction(BigInteger.valueOf(num), BigInteger.valueOf(den)); }
    public static Fraction of(BigInteger k) { return new Fraction(k); }

    /** Parse "a/b" or "a" (whitespace ok). */
    public static Fraction parse(String s) {
        String t = s.trim();
        int slash = t.indexOf('/');
        if (slash < 0) return new Fraction(new BigInteger(t), BigInteger.ONE);
        BigInteger a = new BigInteger(t.substring(0, slash).trim());
        BigInteger b = new BigInteger(t.substring(slash + 1).trim());
        return new Fraction(a, b);
    }

    public BigInteger numerator()   { return n; }
    public BigInteger denominator() { return d; }

    // ---- Numeric ----

    @Override public Fraction add(Fraction o) {
        // (n/d) + (x/y) with cross-cancel to limit growth
        BigInteger g1 = n.gcd(o.d);
        BigInteger g2 = d.gcd(o.n);
        BigInteger a = n.divide(g1);
        BigInteger b = o.n.divide(g2);
        BigInteger c = d.divide(g2);
        BigInteger e = o.d.divide(g1);
        return new Fraction(a.multiply(e).add(b.multiply(c)), c.multiply(e));
    }

    @Override public Fraction subtract(Fraction o) {
        // reuse add on (-o)
        return add(o.negate());
    }

    @Override public Fraction multiply(Fraction o) {
        // (n/d)*(x/y) with cross-cancel
        BigInteger g1 = n.gcd(o.d);
        BigInteger g2 = d.gcd(o.n);
        BigInteger a = n.divide(g1);
        BigInteger b = o.n.divide(g2);
        BigInteger c = d.divide(g2);
        BigInteger e = o.d.divide(g1);
        return new Fraction(a.multiply(b), c.multiply(e));
    }

    @Override public Fraction divide(Fraction o) {
        if (o.n.signum() == 0) throw new ArithmeticException("Divide by zero fraction");
        // (n/d) / (x/y) = (n*y)/(d*x) with cross-cancel
        BigInteger g1 = n.gcd(o.n);
        BigInteger g2 = d.gcd(o.d);
        BigInteger a = n.divide(g1);
        BigInteger b = o.d.divide(g2);
        BigInteger c = d.divide(g2);
        BigInteger e = o.n.divide(g1);
        return new Fraction(a.multiply(b), c.multiply(e));
    }

    public Fraction inverse() {
        if (n.signum() == 0) throw new ArithmeticException("Zero has no inverse");
        return new Fraction(d, n);
    }

    @Override public Fraction negate() { return n.signum() == 0 ? ZERO : new Fraction(n.negate(), d); }
    @Override public Fraction abs()    { return n.signum() < 0 ? negate() : this; }
    @Override public int signum()      { return n.signum(); }
    @Override public boolean isZero()  { return n.signum() == 0; }

    // ---- Comparable ----
    @Override public int compareTo(Fraction o) {
        // a/b ? c/d  <=>  ad ? cb
        return n.multiply(o.d).compareTo(o.n.multiply(d));
    }

    // ---- Object ----
    @Override public boolean equals(Object obj) {
        if (this == obj) return true;
        if (!(obj instanceof Fraction)) return false;
        Fraction o = (Fraction) obj;
        return n.equals(o.n) && d.equals(o.d);
    }

    @Override public int hashCode() { return n.hashCode() * 31 + d.hashCode(); }

    /** lrs-style: integers print as integers, rationals as a/b */
    @Override public String toString() {
        return d.equals(BigInteger.ONE) ? n.toString() : n + "/" + d;
    }
}
