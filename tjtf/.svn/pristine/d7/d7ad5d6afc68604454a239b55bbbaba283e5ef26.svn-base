package javaslam.util;

/**
 * A class containing methods for counting floating point operations.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.3 $ ($Date: 2003/01/07 01:16:46 $)
 */
public class Flops {

  /**
   * A flag that determines whether {@link #count(long)} will
   * increment the flop count.
   */
  public static boolean countFlops = true;

  /**
   * A counter of the number of floating-point operations performed.
   */
  protected static long flops = 0L;

  /**
   * Returns the number of floating point operations performed by code
   * in the <code>tjtf</code> package since the virtual machine was
   * started.
   */
  public static long count() {
    return flops;
  }

  /**
   * Increments the flop counter by <code>f</code>.
   */
  public static void count(long f) {
    if (countFlops)
      flops += f;
  }

  /**
   * A method that resets the flop count.
   */
  public static void reset() {
    flops = 0L;
  }

  /**
   * Counts the number of floating point operations used to solve for
   * <i>X</i> in the matrix equation <i>AX</i> = <i>B</i> where
   * <i>A</i> is possibly symmetric positive definite.
   *
   * This algorithm comes from the Numerical Recipes algorithm via the
   * Lightspeed Matlab library of Tom Minka.  
   *
   * @param n     the number of rows in <i>A</i>
   * @param k     the number of columns in <i>A</i> (and also the number
   *              of rows in <i>B</i>)
   * @param m     the number of columns in <i>B</i>
   * @param isSPD <tt>true</tt> if <i>A</i> is positive definite, enabling 
   *              the Cholesky decomposition
   * @return the exact number of floating point operations required to 
   *         solve for <i>X</i>.
   */
  public static int solve(int n, int k, int m, boolean isSPD) {
    if (n == k) { // square A
      if (n == 1) // scalar division
	return m;
      else if (isSPD)
	return chol(n) + backSubst(n, k, m) + forwardSubst(n, k, m);
      else
	return lu(n) + backSubst(n, k, m) + forwardSubst(n, k, m);
    } else if (n > k)
      // This comes from X = (A'*A)\(A'*B)
      return mult(k, n, k) + mult(k, n, m) + 
	solve(k, k, m, isSPD);
    else
      // This comes from X = A'*(A*A')\B
      return mult(n, k, n) + mult(k, n, m) + 
	solve(n, n, m, isSPD);
  }

  /**
   * Counts the number of floating point operations used to solve for
   * <i>X</i> via back substitution in the matrix equation <i>AX</i> =
   * <i>B</i> where <i>A</i> is upper triangular.
   *
   * This algorithm comes from the Numerical Recipes algorithm via the
   * Lightspeed Matlab library of Tom Minka.  
   *
   * @param n     the number of rows in <i>A</i>
   * @param k     the number of columns in <i>A</i> (and also the number
   *              of rows in <i>B</i>)
   * @param m     the number of columns in <i>B</i>
   * @return the exact number of floating point operations required to 
   *         solve for <i>X</i> via back substitution.  
   */
  public static int backSubst(int n, int k, int m) {
    return n * k * m;
  }

  /**
   * Counts the number of floating point operations used to solve for
   * <i>X</i> via back substitution in the matrix equation <i>AX</i> =
   * <i>B</i> where <i>A</i> is unit lower triangular.
   *
   * This algorithm comes from the Numerical Recipes algorithm via the
   * Lightspeed Matlab library of Tom Minka.  
   *
   * @param n     the number of rows in <i>A</i>
   * @param k     the number of columns in <i>A</i> (and also the number
   *              of rows in <i>B</i>)
   * @param m     the number of columns in <i>B</i>
   * @return the exact number of floating point operations required to 
   *         solve for <i>X</i> via back substitution.  
   */
  public static int forwardSubst(int n, int k, int m) {
    return n * k * m;
  }

  /**
   * Counts the number of floating point operations used to invert an
   * <code>n</code>-by-<code>n</code> matrix that is possibly
   * symmetric positive definite.
   *
   * This algorithm comes from the Numerical Recipes algorithm via the
   * Lightspeed Matlab library of Tom Minka.  
   *
   * @param n     the number of rows and columns of the square matrix
   * @param isSPD <tt>true</tt> if the matrix is positive definite, 
   *              enabling the Cholesky decomposition
   * @return the exact number of floating point operations required to 
   *         invert the matrix
   */
  public static int inv(int n, boolean isSPD) {
    return solve(n, n, n, isSPD);
  }
    
  /**
   * Counts the number of floating point operations used to compute
   * the Cholesky decomposition of an <code>n</code>-by-<code>n</code>
   * symmetric positive definite matrix.
   *
   * This formula comes from the Numerical Recipes algorithm via the
   * Lightspeed Matlab library of Tom Minka.  
   *
   * @param n     the number of rows and columns of the square matrix
   * @return the number of flops required to compute the Cholesky
   *         decomposition 
   */
  public static int chol(int n) {
    return ((2 * n * n * n + 3 * n * n + n) / 6) + n * (Flops.sqrt() - 1);
  }

  /**
   * Counts the number of floating point operations used to compute
   * the LU decomposition of an <code>n</code>-by-<code>n</code>
   * square matrix.
   *
   * This formula comes from 
   *      http://www.maths.uq.edu.au/~gac/math2200/mn_nla2.pdf.
   * 
   * @param n     the number of rows and columns of the square matrix
   * @return the number of flops required to compute the Cholesky
   *         decomposition 
   */
  public static int lu(int n) {
    int k = 2 * (n - 1) * n;
    return ((n * k) / 3) + (k / 12);
  }

  /**
   * Counts the number of floating point operations used to compute
   * the square root of a number.  This function returns
   * <code>15</code> and is based on source code for
   * <code>sqrt()</code> function at
   * http://www.opencores.org/cvsweb.shtml/or1k/newlib/newlib/libm/mathfp/s_sqrt.c.
   * This formula comes from the Lightspeed Matlab library of Tom
   * Minka.
   */
  public static int sqrt() {
    return 15;
  }

  /**
   * Counts the number of floating point operations used to compute
   * the determinant of an <code>n</code>-by-<code>n</code> matrix
   * that is possibly symmetric positive definite.
   *
   * @param n     the number of rows and columns of the square matrix
   * @param isSPD <tt>true</tt> if the matrix is positive definite, 
   *              enabling the Cholesky decomposition
   *
   * This formula comes from the Lightspeed Matlab library of Tom
   * Minka.  
   */
  public static int det(int n, boolean isSPD) {
    if (n == 1) return 1;
    else return n + (isSPD? chol(n) : lu(n));
  }

  /**
   * Counts the number of floating point operations used to compute
   * the logarithm of a number.  This function returns
   * <code>20</code>; this formula comes from the Lightspeed Matlab
   * library of Tom Minka.
   */
  public static int log() {
    return 20;
  }

  /**
   * Counts the number of floating point operations used to compute
   * the exponential of a number.  This function returns
   * <code>20</code>; this formula comes from the Lightspeed Matlab
   * library of Tom Minka.
   */
  public static int exp() {
    return 20;
  }

  /**
   * Counts the number of floating point operations used to compute a
   * random double-precision floating point number between 0.0 and 1.0
   * using a simple linear-congruential formula <i>x</i><sub>i +
   * 1</sub> = <b>a</b> <i>x</i><sub>i</sub> + <b>b</b> (mod <b>c</b>).
   */
  public static int rand() {
    return 3;
  }

  /**
   * Counts the number of floating point operations used to compute a
   * sample from the standard normal using the Box-Muller method.
   */
  public static int randn() {
    // log is similar to cos/sin
    return 2 * rand() + 2 * log() + sqrt() + 1;
  }

  /**
   * Counts the number of floating point operations used to compute
   * <code>n</code> samples from a <code>k</code>-dimensional
   * multivariate Gaussian distribution represented using a covariance
   * matrix.
   *
   * @param n the number of samples
   * @param k the dimension of the Gaussian distribution
   * @return the number of floating-point operations required to
   *         sample <code>n</code> times from a
   *         <code>k</code>-dimensional multivariate Gaussian
   *         distribution represented using a covariance matrix.  
   */
  public static int randnorm(int n, int k) {
    // Sample k * n times from randn, multiply by the Cholesky, and
    // then add in the mean.
    return randn() * k * n + chol(k) + mult(k, k, n) + k * n; 
  }

  /**
   * Counts the number of floating point operations used to compute
   * the product of two matrices <I>A</i> and <i>B</i>.
   *
   * This formula comes from the Lightspeed Matlab library of Tom
   * Minka.  
   *
   * @param n     the number of rows in <i>A</i>
   * @param k     the number of columns in <i>A</i> (and also the number
   *              of rows in <i>B</i>)
   * @param m     the number of columns in <i>B</i>
   * @return the number of flops used to compute the product <i>AB</i>
   */
  public static int mult(int n, int k, int m) {
    return m * n * (2 * k - 1);
  }

  /**
   * Counts the number of floating point operations used to compute
   * the trace of a square matrix <I>A</i>.
   *
   * @param n     the number of rows (and columns) in <i>A</i>
   * @return the number of flops used to compute tr(<i>A</i>)
   */
  public static int trace(int n) {
    return n - 1;
  }

  /**
   * Counts the number of floating point operations used to add two
   * <code>n</code>-by-<code>m</code> matrices.
   *
   * @param n     the number of rows in each matrix
   * @param m     the number of columns in each matrix
   * @return      the number of flops used to compute the sum 
   *              <i>A</i> + <i>B</i>
   */
  public static int add(int n, int m) {
    return m * n;
  }
}
