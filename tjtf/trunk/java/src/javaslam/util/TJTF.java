package javaslam.util;

import Jama.*;

/**
 * A class containing package-level utilities.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.8 $ ($Date: 2003/02/07 20:06:38 $)
 */
public class TJTF {

  /**
   * Determines whether routine status messages should be displayed.
   *
   * @param on the new value
   */
  public static void setVerbose(boolean on) {
    verbose = on;
  }
    
  /**
   * Determines whether mundane status messages should be displayed.
   *
   * @param on the new value
   */
  public static void setBlather(boolean on) {
    blather = on;
  }
    
  /**
   * Determines whether special (expensive) checks should be made to
   * ensure that everything is operating correctly.
   *
   * @param on the new value
   */
  public static void setDebug(boolean on) {
    debug = on;
  }
    
  /**
   * A flag indicating whether routine status messages should be
   * displayed.
   */
  public static boolean verbose = false;

  /**
   * A flag indicating whether mundane status messages should be
   * displayed.
   */
  public static boolean blather = false;

  /**
   * A flag indicating whether special (expensive) checks should be
   * made to ensure that everything is operating correctly.
   */
  public static boolean debug = false;

  /**
   * A method that reports a message.
   *
   * @param s the message to be reported
   */
  public static void say(String s) { System.out.println(s); }

  /**
   * Creates a Matlab string representation of a matrix.
   *
   * @param m the matrix
   * @return a string representing <code>m</code> in Matlab format
   */
  public static String toString(Matrix m) {
    StringBuffer b = new StringBuffer();
    java.text.DecimalFormat f = new java.text.DecimalFormat("0.##E0");
    b.append("[");
    for (int i = 0; i < m.getRowDimension(); i++) {
      for (int j = 0; j < m.getColumnDimension(); j++)
        b.append("\t" + f.format(m.get(i, j)));
      if (i < m.getRowDimension() - 1) 
	b.append("; ...\n");
    }
    b.append("]\n");
    return b.toString();
  }

  /**
   * Transforms an asymmetric matrix <i>M</i> into a symmetric matrix
   * by setting it to (<i>M</i> + <i>M</i><sup>T</sup>)/2.  If
   * <i>M</i> is positive definite, then this method will not alter
   * its eigensystem.
   */
  public static Matrix symmetrize(Matrix m) {
    if (m.getRowDimension() != m.getColumnDimension())
      throw new IllegalArgumentException("Matrix not square");
    for (int i = 0; i < m.getRowDimension(); i++)
      for (int j = i + 1; j < m.getRowDimension(); j++) {
	double d = (m.get(i, j) + m.get(j, i)) / 2.0;
	m.set(i, j, d); 
	m.set(j, i, d); 
      }
    return m;
  }

  /**
   * Regularizes a covariance matrix so that its condition number is
   * no greater than a specified value.  This is accomplished by
   * adding in a small multiple of the identity matrix.
   *
   * @param C    the covariance matrix
   * @param cond the maximum tolerable condition number */
  public static void regularize(Matrix C, double cond) {
    int k = C.getRowDimension();
    // Symmetrize C so its eigenvalues are all real.
    TJTF.symmetrize(C);
    // Compute its eigensystem.
    EigenvalueDecomposition eig = C.eig();
    Matrix D = eig.getD();
    Matrix V = eig.getV();
    // Compute the maximum eigenvalue.
    double maxEig = Double.NEGATIVE_INFINITY;
    double minEig = Double.POSITIVE_INFINITY;
    for (int i = 0; i < k; i++) {
      if (D.get(i, i) > maxEig) maxEig = D.get(i, i);
      if (D.get(i, i) < minEig) minEig = D.get(i, i);
    }
    if (minEig < 0.0d) minEig = 0.0d;
    if (maxEig / minEig > cond) {
      // Compute the minimum permissible eigenvalue.
      minEig = maxEig / cond;
      // Regularize the eigensystem by adding minEig * I.
      for (int i = 0; i < k; i++) 
	D.set(i, i, D.get(i, i) + minEig);
      C.setMatrix(0, k - 1, 0, k - 1, V.times(D.times(V.transpose())));
    }
  }
}
