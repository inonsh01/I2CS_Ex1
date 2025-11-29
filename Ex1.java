//package assignments.Ex1;

import java.util.Arrays;

/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions -
 * represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function:
 * 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe
 * 
 */
public class Ex1 {
  /**
   * Epsilon value for numerical computation, it serves as a "close enough"
   * threshold.
   */
  public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
  /**
   * The zero polynomial function is represented as an array with a single (0)
   * entry.
   */
  public static final double[] ZERO = { 0 };

  /**
   * Computes the f(x) value of the polynomial function at x.
   * 
   * @param poly - polynomial function
   * @param x
   * @return f(x) - the polynomial function value at x.
   */
  public static double f(double[] poly, double x) {
    double ans = 0;
    for (int i = 0; i < poly.length; i++) {
      double c = Math.pow(x, i);
      ans += c * poly[i];
    }
    return ans;
  }

  /**
   * Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
   * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps,
   * assuming p(x1)*p(x2) <= 0.
   * This function should be implemented recursively.
   * 
   * @param p   - the polynomial function
   * @param x1  - minimal value of the range
   * @param x2  - maximal value of the range
   * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
   * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
   */
  public static double root_rec(double[] p, double x1, double x2, double eps) {
    double f1 = f(p, x1);
    double x12 = (x1 + x2) / 2;
    double f12 = f(p, x12);
    if (Math.abs(f12) < eps) {
      return x12;
    }
    if (f12 * f1 <= 0) {
      return root_rec(p, x1, x12, eps);
    } else {
      return root_rec(p, x12, x2, eps);
    }
  }

  /**
   * This function computes a polynomial representation from a set of 2D points on
   * the polynom.
   * The solution is based on: //
   * http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
   * Note: this function only works for a set of points containing up to 3 points,
   * else returns null.
   * 
   * @param xx
   * @param yy
   * @return an array of doubles representing the coefficients of the polynom.
   */
  public static double[] PolynomFromPoints(double[] xx, double[] yy) {
    double[] ans = null;
    int lx = xx.length;
    int ly = yy.length;
    if (xx != null && yy != null && lx == ly && lx > 1 && lx < 4) {
      // add you code below
      // ----------------

      // define every given point in the right place (xx[0] -> x1 etc.)
      double x1 = xx[0], x2 = xx[1], x3 = xx[2], y1 = yy[0], y2 = yy[1], y3 = yy[2];

      // denominator of the equations
      double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);

      // if the denominator is too close to 0 return null
      if (Math.abs(denom) < EPS) {
        return null;
      }

      // algorithem to calculate A,B,C.
      // to be honest i understood the code but not really the math behind it.
      double A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
      double B = (Math.pow(x3, 2) * (y1 - y2) + Math.pow(x2, 2) * (y3 - y1) + Math.pow(x1, 2) * (y2 - y3)) / denom;
      double C = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

      ans = new double[] { C, B, A };

      // to avoid cases like -0.0 etc.
      for (int i = 0; i < ans.length; i++) {
        if (Math.abs(ans[i]) < EPS) {
          ans[i] = 0;
        }
      }

      // ----------------
    }

    return ans;
  }

  /**
   * Two polynomials functions are equal if and only if they have the same values
   * f(x) for n+1 values of x,
   * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
   * 
   * @param p1 first polynomial function
   * @param p2 second polynomial function
   * @return true iff p1 represents the same polynomial function as p2.
   */
  public static boolean equals(double[] p1, double[] p2) {
    boolean ans = true;
    // add you code below
    // -------------------

    // checking edge cases

    // if they both empty they equals by definition
    if (p1 == null && p2 == null) {
      return true; // they are equals, return true.
    }

    // get the shorter polynom
    int minLength = Math.min(p1.length, p2.length);
    // loop on every index until the end of the shorter and compare them, if one of
    // them not equals to other -
    // they are not equals
    for (int i = 0; i < minLength; i++) {
      // if only one of them (xor) zero they are not in the same degree so they are
      // not equals.
      // or if they are in the same degree check that the diffrence not bigger then
      // epsilon.
      if ((p1[i] == 0 && p2[i] != 0) || (p1[i] != 0 && p2[i] == 0) || (Math.abs(p1[i] - p2[i]) > EPS)) {
        ans = false;
      }
    }

    // if the length of each not equals to the other
    if (p1.length != p2.length) {

      // if the rest of the longer polynom is zeroes they are still equals -
      // for example: p([1,2]) = p([1,2,0,0]).
      // so we need to check if all the indexes are zeroes in the longer polynom
      // after the length of the short one.
      if (!isAllZero((p1.length > p2.length) ? p1 : p2, minLength)) {
        return false;
      }
    }

    // ------------------
    return ans;
  }

  /**
   * Computes a String representing the polynomial function.
   * For example the array {2,0,3.1,-1.2} will be presented as the following
   * String "-1.2x^3 +3.1x^2 +2.0"
   * 
   * @param poly the polynomial function represented as an array of doubles
   * @return String representing the polynomial function:
   */
  public static String poly(double[] poly) {
    String ans = "";
    if (poly.length == 0) {
      ans = "0";
    } else {
      // add you code below
      // -------------------

      // loop over the polynom from the last one to the first
      for (int i = poly.length - 1; i >= 0; i--) {
        if (i == 0) {
          // the last (x ^ 0) will be just the number because it multiple by 1
          ans += poly[i];
        } else if (i == 1) {
          // the (x ^ 1) will be just the number multuple with x
          ans += poly[i] + "x ";
        } else {
          // other will be with multiple by x^i
          ans += poly[i] + "x^" + i + " ";
        }
        // if the next is positive add "+" (otherwise keep the minus)
        if (i > 1 && poly[i - 1] > 0) {
          ans += "+";
        }

      }

      // -------------------
    }
    return ans;
  }

  /**
   * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps.
   * This function computes an x value (x1<=x<=x2)
   * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <=
   * 0.
   * 
   * @param p1  - first polynomial function
   * @param p2  - second polynomial function
   * @param x1  - minimal value of the range
   * @param x2  - maximal value of the range
   * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
   * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
   */
  public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
    double ans = x1;
    // add you code below
    // ----------------

    // To find where two polynomails is with the same value
    // (up to epsilon) at some range,
    // we can subtract them and check when the value of
    // subtracted polynomail at x is close to epsilon (root).
    double[] subPol = subtract(p1, p2);
    ans = root_rec(subPol, x1, x2, eps);

    // ----------------
    return ans;
  }

  /**
   * Given a polynomial function (p), a range [x1,x2] and an integer with the
   * number (n) of sample points.
   * This function computes an approximation of the length of the function between
   * f(x1) and f(x2)
   * using n inner sample points and computing the segment-path between them.
   * assuming x1 < x2.
   * This function should be implemented iteratively (none recursive).
   * 
   * @param p                - the polynomial function
   * @param x1               - minimal value of the range
   * @param x2               - maximal value of the range
   * @param numberOfSegments - (A positive integer value (1,2,...).
   * @return the length approximation of the function between f(x1) and f(x2).
   */
  public static double length(double[] p, double x1, double x2, int numberOfSegments) {
    double ans = x1;
    // add you code below
    // -----------------

    double x;
    // reset ans
    ans = 0.0;

    // define step by length of the range devide by number of segments
    double step = Math.abs(x2 - x1) / numberOfSegments;

    // loop over every segment
    for (int i = 0; i < numberOfSegments; i++) {
      // current x is the start point + number of steps
      // depends on current index of segment
      x = x1 + (step * i);
      // calculate length (using pythagoras) of segment and add it to answer
      ans += Math.sqrt(Math.pow(step, 2) + Math.pow(f(p, x + step) - f(p, x), 2));
    }

    // -----------------
    return ans;
  }

  /**
   * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer
   * representing the number of Trapezoids between the functions (number of
   * samples in on each polynom).
   * This function computes an approximation of the area between the polynomial
   * functions within the x-range.
   * The area is computed using Riemann's like integral
   * (https://en.wikipedia.org/wiki/Riemann_integral)
   * 
   * @param p1                - first polynomial function
   * @param p2                - second polynomial function
   * @param x1                - minimal value of the range
   * @param x2                - maximal value of the range
   * @param numberOfTrapezoid - a natural number representing the number of
   *                          Trapezoids between x1 and x2.
   * @return the approximated area between the two polynomial functions within the
   *         [x1,x2] range.
   */
  public static double area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid) {
    double ans = 0;
    // add you code below
    // ---------------

    double currX, nextX;
    // define step by length of the range devide by number of Trapezoids
    double step = Math.abs(x2 - x1) / numberOfTrapezoid;

    // find the diffrence polynomial function
    double[] subPol = subtract(p1, p2);

    // loop over every segment
    for (int i = 0; i < numberOfTrapezoid; i++) {
      // define this x and the next ( + step)
      currX = x1 + (step * i);
      nextX = currX + step;

      // if there's a sign change (root crossing)
      if (f(subPol, currX) * f(subPol, nextX) <= 0) {

        // find the root of the substract polynomial function
        double r = root_rec(subPol, currX, nextX, 0.0001);

        // calculate from current x till root
        ans += generalArea(p1, p2, currX, r);
        // calculate from root till x + step
        ans += generalArea(p1, p2, r, nextX);

        // if current segment not crossing root, calculate as usual
      } else {
        ans += generalArea(p1, p2, currX, nextX);
      }
    }

    // ---------------
    return ans;
  }

  /**
   * This function computes the array representation of a polynomial function from
   * a String
   * representation. Note:given a polynomial function represented as a double
   * array,
   * getPolynomFromString(poly(p)) should return an array equals to p.
   * 
   * @param p - a String representing polynomial function.
   * @return
   */
  public static double[] getPolynomFromString(String p) {
    double[] ans = ZERO;// -1.0x^2 +3.0x +2.0
    // add you code below
    // ------------------

    // split the sring by spaces to create array of values of the polynom
    String[] strArr = p.split("\\s+");

    // assuming that the degrees writing by order from the biggest to the smallest
    // get the max degree by extracting the first element degree
    int maxDegree = Integer.parseInt(strArr[0].split("x\\^")[1]);

    // create new array with the size of to max degree + 1 (because the 0 degree)
    ans = new double[maxDegree + 1];

    // create integer to keep the current index of the answer array.
    int ansCurrIndex = 0;

    // loop all over the polynom in reverse order
    for (int i = strArr.length - 1; i >= 0; i--) {

      // split every element by x and ^ to split the coefficient of x and the degree
      String[] polItem = strArr[i].split("x\\^");

      // if degree is 0 or 1
      if (polItem.length == 1) {
        // if degree is 1
        if (polItem[0].contains("x")) {
          // if 1 degree is at index 0 add to the first 0
          if (ansCurrIndex == 0) {
            ans[0] = 0;
          }
          // insert the coefficient without "x" or "+" or spaces
          ans[1] = Double.parseDouble(polItem[0].replaceAll("[ +x]", ""));
          // promote two steps forward
          ansCurrIndex = 2;

          // if degree is 0
        } else {
          // insert just the variable
          ans[0] = Double.parseDouble(polItem[0]);
          // promote one step forward
          ansCurrIndex = 1;
        }

        // if degree is bigger then 1
      } else {
        // casting for current degree and current coefficient (without "+" or spaces)
        double currCeo = Double.parseDouble(polItem[0].replaceAll("[ +]", ""));
        int currDeg = Integer.parseInt(polItem[1]);

        // if the degree is in the right index just insert the coefficient
        if (ansCurrIndex == currDeg) {
          ans[ansCurrIndex] = currCeo;
        } else {
          // loop from answer current index until current degree
          for (int j = ansCurrIndex; j < currDeg + 1; j++) {
            // if currDeg is at the right index insert it
            if (j == currDeg) {
              ans[ansCurrIndex] = currCeo;

              // otherwise put zero in every cell
            } else {
              ans[ansCurrIndex] = 0;
            }

            // promote one step forward
            ansCurrIndex++;
          }
        }
      }
    }

    // ------------------
    return ans;
  }

  /**
   * This function computes the polynomial function which is the sum of two
   * polynomial functions (p1,p2)
   * 
   * @param p1 - first polynomial function
   * @param p2 - second polynomial function
   * @return - the sum function of the polynomails
   */
  public static double[] add(double[] p1, double[] p2) {
    double[] ans = ZERO;
    // add you code below
    // ---------------

    // define the length of the longer polynom
    int maxLength = Math.max(p1.length, p2.length);

    // create new array with the correct size
    ans = new double[maxLength];

    // loop over the length of the longer polynom
    for (int i = 0; i < maxLength; i++) {
      // if p1[i] not exist insert just p2[i] to ans[i]
      if (p1.length - 1 < i) {
        ans[i] = p2[i];

        // if p2[i] not exist insert just p1[i] to ans[i]
      } else if (p2.length - 1 < i) {
        ans[i] = p1[i];

        // if p1[i] and p2[i] exist just add them
      } else {
        ans[i] = p1[i] + p2[i];
      }
    }

    // ---------------
    return ans;
  }

  /**
   * This function computes the polynomial function which is the multiplication of
   * two polynoms (p1,p2)
   * 
   * @param p1 - first polynomial function
   * @param p2 - second polynomial function
   * @return - the multiplication function
   */
  public static double[] mul(double[] p1, double[] p2) {
    double[] ans = ZERO;//
    // add you code below
    // ---------------

    // create new array with the correct size
    ans = new double[p1.length + p2.length - 1];

    // loop over all polynom 1
    for (int i = 0; i < p1.length; i++) {
      // loop over all polynom 2
      for (int j = 0; j < p2.length; j++) {
        // multiply p[i] and p2[j] and add it to corrent ans[i+j]
        ans[i + j] += p1[i] * p2[j];
      }
    }

    // ---------------
    return ans;
  }

  /**
   * This function computes the derivative of the p0 polynomial function.
   * 
   * @param po - polynomial function
   * @return - the derivative function
   */
  public static double[] derivative(double[] po) {
    double[] ans = ZERO;//
    // add you code below
    // ---------------

    // edge cases:
    // assuming that po is not empty
    // otherwise if(po.length == 0) return [0]
    // if the polynom's max degree is 0 return just [0]
    if (po.length == 1) {
      return ans;
    }

    // otherwise create new array with correct length (while derivativing
    // the length of degrees drops by 1)
    ans = new double[po.length - 1];

    // skip above the first one (at 0 degrees)
    for (int i = 1; i < po.length; i++) {
      // insert to every cell the derivative of the next index. for example
      // [a,x,2x,2x] -> [x,4x,6x]
      ans[i - 1] = po[i] * i;
    }

    // --------------
    return ans;
  }

  /**
   * Check if all the values in the array are zero
   * 
   * @param po the polynomial function represented as an array of doubles
   * @param in the index to start checking from
   * @return Boolean value indicating whether all values in the array are zero or
   *         not
   */
  public static boolean isAllZero(double[] po, int in) {
    // looping over the array from give index
    for (int i = in; i < po.length; i++) {
      // if exist value that isn't 0 return false
      if (po[i] != 0) {
        return false;
      }
    }
    // all the values are 0, return true
    return true;
  }

  /**
   * * Given two polynomial functions (p1,p2), a range [x1,x2]
   * This function is general function to computes an approximation of the area
   * between the polynomial
   * functions within the x-range.
   * 
   * @param p1 - first polynomial function
   * @param p2 - second polynomial function
   * @param x1 - start value of x to calculate
   * @param x2 - start value of xto calculate
   * @return the approximated area between the two polynomial functions within the
   *         [x1,x2] range.
   */
  public static double generalArea(double[] p1, double[] p2, double x1, double x2) {
    double ans = 0;

    // calc width - distance between x-ses
    double width = Math.abs(x2 - x1);

    // find first point's y minus the second at absolute value it gives the height
    double h1 = Math.abs(f(p1, x1) - f(p2, x1));
    double h2 = Math.abs(f(p1, x2) - f(p2, x2));

    // algorithem to find area of Trapezoid is (height's sum * width / 2)
    // here the step defind as a width
    ans += (h1 + h2) * width / 2.0;
    return ans;
  }

  /**
   * This function computes the polynomial function which is the subtract of two
   * polynomial functions (p1,p2)
   * 
   * @param p1 - first polynomial function
   * @param p2 - second polynomial function
   * @return the subtract function of the polynomails
   */
  public static double[] subtract(double[] p1, double[] p2) {

    double ceo1, ceo2;
    // define the lenght by the longer polynom
    int length = Math.max(p1.length, p2.length);

    // define new answer array with max length
    double[] ans = new double[length];

    // loop over the polynoms
    for (int i = 0; i < length; i++) {

      // inline way to insert into coefficient -
      // if index smaller them the max degree, return coefficient in this index.
      // else return 0.0 (0*^i)
      ceo1 = (i < p1.length) ? p1[i] : 0.0;
      ceo2 = (i < p2.length) ? p2[i] : 0.0;

      // subtract coefficients
      ans[i] = ceo1 - ceo2;
    }
    return ans;
  }
}