package javaslam.filter;

import Jama.*;
import javaslam.prob.*;

/**
 * A vector-valued function that takes a vector input, the bottom part
 * of which is a white noise vector.  In particular, if
 * <code>dim</code> is the result of calling {@link
 * Gaussian#getDimension()} on the Gaussian {@link #getNoiseModel()},
 * then elements <code>{@link VectorFunction#getInputDim()} -
 * dim</code> through <code>{@link VectorFunction#getInputDim()} -
 * 1</code> of the input vector are the noise variables.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.1 $ ($Date: 2003/02/07 02:25:45 $) */
public class LinearGaussianFunction 
  implements NoisyVectorFunction, ExtendedVectorFunction {

  /**
   * The linear coefficient.
   */
  protected Matrix A;

  /**
   * The noise model.
   */
  protected Gaussian q;

  /**
   * Constructor.
   *
   * @param A the linear coefficient
   * @param q a Gaussian noise model (in moment form) that has
   *          as many dimensions as A has columns
   */
  public LinearGaussianFunction(Matrix A, Gaussian q) {
    this.A = A;
    this.q = q;
    if (A.getRowDimension() != q.getDimension())
      throw new IllegalArgumentException("The dimension of the noise " + 
					 "model must match the number of " + 
					 "rows of the matrix coefficient");
  }

  /**
   * Returns the input dimension of this function.
   */
  public int getInputDim() {
    return 2 * A.getColumnDimension();
  }

  /**
   * Returns the output dimension of this function.
   */
  public int getOutputDim()  {
    return A.getRowDimension();
  }

  /**
   * Evaluates this function at the supplied input.
   *
   * @param input an array with {@link #getInputDim()} elements
   * @return an array with {@link #getOutputDim()} elements
   */
  public double[] evaluate(double[] input) {
    Matrix tmp = new Matrix(input, input.length);
    Matrix x = tmp.getMatrix(0, A.getColumnDimension() - 1, 0, 0);
    Matrix w = tmp.getMatrix(A.getColumnDimension(), getInputDim() - 1, 0, 0);
    Matrix out = A.times(x).plus(w);
    return out.getColumnPackedCopy();
  }

  /**
   * Returns the Gaussian distribution over the noise input.
   *
   * @return the Gaussian distribution over the noise input
   */
  public Gaussian getNoiseModel() {
    return new Gaussian(q);
  }

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
  public double[][] jacobian(double[] input) {
    int d = A.getColumnDimension();
    Matrix J = new Matrix(getOutputDim(), getInputDim());
    J.setMatrix(0, getOutputDim() - 1, 0, d - 1, A);
    J.setMatrix(0, getOutputDim() - 1, d, getInputDim() - 1, 
		Matrix.identity(d, d));
    return J.getArray();
  }
}
