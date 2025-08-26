# Vertex Enumeration in Java

This repository provides a clean-room Java implementation of the
vertex‑enumeration algorithm originally described by David Avis and
Komei Fukuda in their 1992 paper
*“A pivoting algorithm for convex hulls and vertex enumeration of
arrangements and polyhedra”*.  It is inspired by the C codebase found in
the `lrslib` project but does not copy or reuse any of that source.
Instead, the code here follows the structure of the original algorithm
and uses modern Java idioms.

## Project structure

```
vertex-enumeration-java/
├── README.md           – This file
├── pom.xml             – Maven build file (Java 17)
└── src/main/java/com/vertexenumeration/
    ├── Fraction.java   – Arbitrary‑precision rational numbers
    ├── Matrix.java     – Simple dense matrix of Fraction values
    ├── Polyhedron.java – Parse/write H‑ and V‑representations
    ├── Dictionary.java – Tableau used by the reverse search
    ├── VertexEnumerator.java – High‑level enumeration algorithm
    └── Main.java       – Command‑line entry point
```

### Fraction

`Fraction` wraps a numerator and denominator stored as `BigInteger` and
reduces every fraction to lowest terms.  Basic arithmetic (addition,
subtraction, multiplication and division) is provided.  Overflow is
avoided by using Java’s arbitrary‑precision integers.

### Matrix

`Matrix` represents a two‑dimensional array of `Fraction` values.  It is
used for storing the constraint matrix of an H‑representation or the
vertex/ray matrix of a V‑representation.  The class provides methods
for reading from and writing to the plain text file format used by
lrslib; see `Polyhedron` for higher level parsing logic.

### Polyhedron

`Polyhedron` encapsulates either an H‑representation or a V‑representation
of a convex polyhedron.  It contains the input rows, the type of
representation, and helper methods to read from the lrslib text format
(with `H‑representation`/`V‑representation` headers) and to write a
representation back out.  This class does not perform any geometry – it
merely holds the data.

### Dictionary

The reverse search algorithm operates on dictionaries (also called
tableaux) that track a basis and cobasis along with the current
solution.  The `Dictionary` class stores the coefficients and keeps
track of which columns are basic versus non‑basic.  Pivoting operations
are provided to walk through adjacent bases according to the
lexicographic pivot rule.  Only the skeleton of the data structure is
implemented here; details of the pivot rule and reverse search must be
filled in when translating the algorithm.

### VertexEnumerator

`VertexEnumerator` ties everything together.  It accepts a
`Polyhedron` in one representation (H or V), constructs an initial
`Dictionary`, and then performs a reverse search to enumerate all
vertices (and rays).  The current implementation contains only method
stubs; the translation of the Avis–Fukuda algorithm should be added
incrementally.  See the paper and the lrslib manual for guidance on
pivot rules, dictionary updates and output handling.

### Main

`Main` is a simple command‑line front‑end.  It reads a file in the
lrslib format, detects whether it is an H‑ or V‑representation, and
invokes `VertexEnumerator` to convert to the dual representation.  It
writes the result to standard output.  Only basic argument parsing is
provided.

## Building and running

This project uses Maven with Java 17.  To compile and run from the
command line:

```sh
mvn package
java -cp target/vertex-enumeration-java-1.0-SNAPSHOT.jar \
    com.vertexenumeration.Main path/to/input.ine
```

Replace `input.ine` with the path to your H‑ or V‑representation file.

## Next steps

1. **Implement the enumeration algorithm.**  The core algorithm is
   described in Avis–Fukuda (1992) and implemented in `lrslib`.  Use
   the provided class stubs to build your own implementation in
   Java.  Pay particular attention to the lexicographic pivot rule
   and the handling of fractional arithmetic.
2. **Add redundancy removal and minimum representations.**  lrslib
   contains routines to detect and remove redundant inequalities and to
   find hidden linearities.  Once vertex enumeration is working,
   consider porting these features.
3. **Create unit tests.**  Sample inputs such as `cube.ine` and
   `cube.ext` can be used to check that the implementation produces
   correct vertices and rays.  A simple JUnit test framework can
   automate these checks.

Contributions and pull requests are welcome.  Please respect the
original `lrslib` license when re‑using ideas or algorithms.