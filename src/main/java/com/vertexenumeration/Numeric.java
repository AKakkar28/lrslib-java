package com.vertexenumeration;

/** Minimal exact-number abstraction for algorithm code. */
public interface Numeric<T extends Numeric<T>> extends Comparable<T> {
    T add(T o);
    T subtract(T o);
    T multiply(T o);
    T divide(T o);
    T negate();
    T abs();
    int signum();
    boolean isZero();
}
