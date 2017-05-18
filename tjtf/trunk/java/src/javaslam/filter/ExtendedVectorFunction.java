package javaslam.filter;

/**
 * A vector-valued function that takes a vector input and which can
 * compute its Jacobian at any input.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.1 $ ($Date: 2003/02/07 02:25:45 $) 
 */
public interface ExtendedVectorFunction extends VectorFunction {

  /**
   * Evaluates the Jacobian of this function at the supplied input.
   *
   * @param input a vector with {@link #getInputDim()} elements
   * @return a matrix with {@link #getOutputDim()} rows and
   *         {@link #getInputDim()} columns; if <b>y</b> = <i>f</i>(<b>x</b>)
   *         then element (<i>i</i>, <i>j</i>) is the partial derivative 
   *         of <b>y</b><sub><i>i</i></sub> with respect to 
   *         <b>x</b><sub><i>j</i></sub> at <code>input</code>
   */
  public double[][] jacobian(double[] input);
}
