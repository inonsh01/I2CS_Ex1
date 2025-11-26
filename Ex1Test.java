import org.junit.jupiter.api.Test;
import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

/**
 * * Introduction to Computer Science 2026, Ariel University,
 * * Ex1: arrays, static functions and JUnit
 *
 * This JUnit class represents a JUnit (unit testing) for Ex1-
 * It contains few testing functions for the polynomial functions as define in
 * Ex1.
 * Note: you should add additional JUnit testing functions to this class.
 *
 * @author boaz.ben-moshe
 */

class Ex1Test {
	static final double[] P1 = { 2, 0, 3, -1, 0 }, P2 = { 0.1, 0, 1, 0.1, 3 };
	static double[] po1 = { 2, 2 }, po2 = { -3, 0.61, 0.2 };
	static double[] po3 = { 2, 1, -0.7, -0.02, 0.02 };
	static double[] po4 = { -3, 0.61, 0.2 };
	static double[] po5 = { -3.93, 0.12, 5, 4 };
	static double[] po6 = { -3.93, 0.12, 5, 4 };
	static double[] po7 = {};
	static double[] po9 = null;

	@Test
	/**
	 * Tests that f(x) == poly(x).
	 */
	void testF() {
		double fx0 = Ex1.f(po1, 0);
		double fx1 = Ex1.f(po1, 1);
		double fx2 = Ex1.f(po1, 2);
		assertEquals(fx0, 2, Ex1.EPS);
		assertEquals(fx1, 4, Ex1.EPS);
		assertEquals(fx2, 6, Ex1.EPS);
	}

	@Test
	/**
	 * Tests that p1(x) + p2(x) == (p1+p2)(x)
	 */
	void testF2() {
		double x = Math.PI;
		double[] po12 = Ex1.add(po1, po2);
		double f1x = Ex1.f(po1, x);
		double f2x = Ex1.f(po2, x);
		double f12x = Ex1.f(po12, x);
		assertEquals(f1x + f2x, f12x, Ex1.EPS);
	}

	@Test
	/**
	 * Tests that p1+p2+ (-1*p2) == p1
	 */
	void testAdd() {
		double[] p12 = Ex1.add(po1, po2);
		double[] minus1 = { -1 };
		double[] pp2 = Ex1.mul(po2, minus1);
		double[] p1 = Ex1.add(p12, pp2);
		assertTrue(Ex1.equals(p1, po1));
	}

	@Test
	/**
	 * Tests that p1+p2 == p2+p1
	 */
	void testAdd2() {
		double[] p12 = Ex1.add(po1, po2);
		double[] p21 = Ex1.add(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}

	@Test
	/**
	 * Tests that p1+0 == p1
	 */
	void testAdd3() {
		double[] p1 = Ex1.add(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, po1));
	}

	@Test
	/**
	 * Tests that p1*0 == 0
	 */
	void testMul1() {
		double[] p1 = Ex1.mul(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, Ex1.ZERO));
	}

	@Test
	/**
	 * Tests that p1*p2 == p2*p1
	 */
	void testMul2() {
		double[] p12 = Ex1.mul(po1, po2);
		double[] p21 = Ex1.mul(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}

	@Test
	/**
	 * Tests that p1(x) * p2(x) = (p1*p2)(x),
	 */
	void testMulDoubleArrayDoubleArray() {
		double[] xx = { 0, 1, 2, 3, 4.1, -15.2222 };
		double[] p12 = Ex1.mul(po1, po2);
		for (int i = 0; i < xx.length; i = i + 1) {
			double x = xx[i];
			double f1x = Ex1.f(po1, x);
			double f2x = Ex1.f(po2, x);
			double f12x = Ex1.f(p12, x);
			assertEquals(f12x, f1x * f2x, Ex1.EPS);
		}
	}

	@Test
	/**
	 * Tests a simple derivative examples - till ZERO.
	 */
	void testDerivativeArrayDoubleArray() {
		double[] p = { 1, 2, 3 }; // 3X^2+2x+1
		double[] pt = { 2, 6 }; // 6x+2
		double[] dp1 = Ex1.derivative(p); // 2x + 6
		double[] dp2 = Ex1.derivative(dp1); // 2
		double[] dp3 = Ex1.derivative(dp2); // 0
		double[] dp4 = Ex1.derivative(dp3); // 0
		assertTrue(Ex1.equals(dp1, pt));
		assertTrue(Ex1.equals(Ex1.ZERO, dp3));
		assertTrue(Ex1.equals(dp4, dp3));
	}

	@Test
	/**
	 * Tests the parsing of a polynom in a String like form.
	 */
	public void testFromString() {
		double[] p = { -1.1, 2.3, 3.1 }; // 3.1X^2+ 2.3x -1.1
		double[] p4 = { -1.1, 0, 0, 2.3, 0, 3.1 };
		String sp2 = "3.1x^2 +2.3x -1.1";

		// add another check - if the polynimal not have the all "middle" degrees -
		// (if the max degree is 5, what happends if ^4 coefficient is zero? etc.)
		String sp3 = "3.1x^5 +2.3x^3 -1.1";
		String sp = Ex1.poly(p);
		double[] p1 = Ex1.getPolynomFromString(sp);
		double[] p2 = Ex1.getPolynomFromString(sp2);
		double[] p3 = Ex1.getPolynomFromString(sp3);
		boolean isSame1 = Ex1.equals(p1, p);
		boolean isSame2 = Ex1.equals(p2, p);
		boolean isSame3 = Ex1.equals(p3, p4);

		if (!isSame1) {
			fail();
		}
		if (!isSame2) {
			fail();
		}
		if (!isSame3) {
			fail();
		}
		// assertEquals(sp, Ex1.poly(p1));
	}

	@Test
	/**
	 * Tests the equality of pairs of arrays.
	 */
	public void testEquals() {
		double[][] d1 = { { 0 }, { 1 }, { 1, 2, 0, 0 } };
		double[][] d2 = { Ex1.ZERO, { 1 + Ex1.EPS / 2 }, { 1, 2 } };
		double[][] xx = { { -2 * Ex1.EPS }, { 1 + Ex1.EPS * 1.2 }, { 1, 2, Ex1.EPS / 2 } };
		for (int i = 0; i < d1.length; i = i + 1) {
			assertTrue(Ex1.equals(d1[i], d2[i]));
		}
		for (int i = 0; i < d1.length; i = i + 1) {
			assertFalse(Ex1.equals(d1[i], xx[i]));
		}
	}

	@Test
	/**
	 * Tests is the sameValue function is symmetric.
	 */
	public void testSameValue2() {
		double x1 = -4, x2 = 0;
		double rs1 = Ex1.sameValue(po1, po2, x1, x2, Ex1.EPS);
		double rs2 = Ex1.sameValue(po2, po1, x1, x2, Ex1.EPS);
		assertEquals(rs1, rs2, Ex1.EPS);
	}

	@Test
	/**
	 * Test the area function - it should be symmetric.
	 */
	public void testArea() {
		double x1 = -4, x2 = 0;
		double a1 = Ex1.area(po1, po2, x1, x2, 100);
		double a2 = Ex1.area(po2, po1, x1, x2, 100);
		assertEquals(a1, a2, Ex1.EPS);
	}

	@Test
	/**
	 * Test the area f1(x)=0, f2(x)=x;
	 */
	public void testArea2() {
		double[] po_a = Ex1.ZERO;
		double[] po_b = { 0, 1 };
		double x1 = -1;
		double x2 = 2;
		double a1 = Ex1.area(po_a, po_b, x1, x2, 1);
		double a2 = Ex1.area(po_a, po_b, x1, x2, 2);
		double a3 = Ex1.area(po_a, po_b, x1, x2, 3);
		double a100 = Ex1.area(po_a, po_b, x1, x2, 100);
		double area = 2.5;
		assertEquals(a1, area, Ex1.EPS);
		assertEquals(a2, area, Ex1.EPS);
		assertEquals(a3, area, Ex1.EPS);
		assertEquals(a100, area, Ex1.EPS);
	}

	@Test
	/**
	 * Test the area function.
	 */
	public void testArea3() {
		double[] po_a = { 2, 1, -0.7, -0.02, 0.02 };
		double[] po_b = { 6, 0.1, -0.2 };
		double x1 = Ex1.sameValue(po_a, po_b, -10, -5, Ex1.EPS);
		double a1 = Ex1.area(po_a, po_b, x1, 6, 8);
		double area = 58.5658;
		assertEquals(a1, area, Ex1.EPS);
	}

	@Test
	/**
	 * Test the isAllZero function.
	 */
	public void testIsAllZero() {
		// define few polynomail function
		double[] po_a = { 0, 0, 0, 0, 0, 0, 0, 0 };
		double[] po_b = { 0, 0, 0 };
		double[] po_c = { 0, 0, 0, 0, 1, 0, 0, 0 };
		double[] po_d = { 0, 0, -1 };

		// call to isAllZero function with the polynomial function
		boolean a1 = Ex1.isAllZero(po_a, 0);
		boolean a2 = Ex1.isAllZero(po_b, 0);
		boolean a3 = Ex1.isAllZero(po_c, 0);
		boolean a4 = Ex1.isAllZero(po_d, 0);

		// test true if everything is zero, otherwise test is false.
		assertTrue(a1);
		assertTrue(a2);
		assertFalse(a3);
		assertFalse(a4);
	}

	@Test
	/**
	 * Test the PolynomFromPoints function.
	 */
	public void testPolynomFromPoints() {
		// define first 3 of 2D points and polynomial function
		double[] xx_a = { 0, 1, -1 };
		double[] yy_a = { 0, 1, 1 };
		double[] pol_a = { 0.0, 0.0, 1.0 };

		// define second 3 of 2D points and polynomial function
		double[] xx_b = { 0, 1, 4 };
		double[] yy_b = { 3, 0, 3 };
		double[] pol_b = { 3, -4, 1 };

		// define thhird 3 of 2D points not in the same length to get null
		double[] xx_c = { 0, 1, -1 };
		double[] yy_c = { 0, 1, 1, 4 };

		// calls to PolynomFromPoints function with the arguments
		double[] p1 = Ex1.PolynomFromPoints(xx_a, yy_a);
		double[] p2 = Ex1.PolynomFromPoints(xx_b, yy_b);
		double[] p3 = Ex1.PolynomFromPoints(xx_c, yy_c);

		// test everyone to it pair
		assertArrayEquals(p1, pol_a, Ex1.EPS);
		assertArrayEquals(p2, pol_b, Ex1.EPS);
		assertEquals(p3, null);
	}

	@Test
	/**
	 * Test the length function.
	 */
	public void length() {

		// define examples returns from length function
		double a1 = Ex1.length(po3, 4, 20, 8);
		double a2 = Ex1.length(po3, 3, 19, 30);
		double a3 = Ex1.length(po3, -4, 2, 40);
		double a4 = Ex1.length(po3, -100, 100, 10);

		// internet calculation
		double r1 = 2783.8097;
		double r2 = 2241.2346;
		double r3 = 12.2471;
		double r4 = 3986000.1488;

		// check each with is pair
		assertEquals(a1, r1, Ex1.EPS);
		assertEquals(a2, r2, Ex1.EPS);
		assertEquals(a3, r3, Ex1.EPS);
		assertEquals(a4, r4, Ex1.EPS);
	}
}