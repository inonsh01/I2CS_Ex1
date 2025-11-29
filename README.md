# I2CS_Ex1

Compute Introduction - Exercise 1

## Overview

This project implements a set of static methods for polynomial function operations. Polynomials are represented as arrays of doubles where the array `{0.1, 0, -3, 0.2}` represents the polynomial function: `0.2x^3 -3x^2 + 0x^1 + 0.1`.

## Files

- `Ex1.java` - Main class with polynomial operations
- `Ex1Test.java` - JUnit tests
- `Ex1_GUI.java` - Graphical user interface
- `StdDraw.java` - Draw

## Functions Implemented

### Core Operations

- `f(poly, x)` - Calculate polynomial at point x
- `add(p1, p2)` - Sum two polynomials
- `mul(p1, p2)` - Multiply two polynomials
- `derivative(p)` - Compute derivative of polynomial

### String Operations

- `poly(p)` - Convert polynomial to string representation
- `getPolynomFromString(s)` - Parse polynomial from string

### Advanced Functions

- `root_rec(p, x1, x2, eps)` - Find root recursivly
- `sameValue(p1, p2, x1, x2, eps)` - Find x where (p1(x) - p2(x)) < eps
- `area(p1, p2, x1, x2, n)` - Calculate area between polynomials
- `length(p, x1, x2, n)` - Approximate length of polynomial

### Utility Functions

- `equals(p1, p2)` - Compare polynomials up to epsilon
- `PolynomFromPoints(xx, yy)` - Get polynomial from points

## Dependencies

- JUnit 5

All JAR files are included in the `lib/` directory.

## Images

- Ex1_GUI output

![EX1 GUI](./GUI.png)
