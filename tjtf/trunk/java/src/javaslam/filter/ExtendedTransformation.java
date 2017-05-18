package javaslam.filter;

import java.util.*;
import Jama.*;
import javaslam.prob.*;
import javaslam.util.*;

/**
 * The extended transformation for linearizing functions with
 * Gaussian-distributed inputs.  This technique computes the
 * first-order Taylor expansion of the nonlinear function.  The
 * function to be linearized must implement {@link
 * ExtendedVectorFunction} so that the Jacobian can be computed.
 *
 * <p>This class records counts of all floating point operations using
 * {@link Flops#count(long)}.</p>
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.1 $ ($Date: 2003/02/07 02:25:45 $) 
 */
public class ExtendedTransformation extends Linearization {

  /**
   * A factory for creating linearizations using the extended
   * transformation.
   */
  protected static LinearizationFactory factory = 
    new LinearizationFactory() {
      public Linearization linearize(NoisyVectorFunction f, Gaussian px) {
	return new ExtendedTransformation(f, px);
      }
    };

  /**
   * Constructor.
   *
   * @param f  an {@link ExtendedVectorFunction} that takes a 
   *           <nobr>(<i>n</i> + <i>m</i>)</nobr>-by-1 vector as input 
   *           and returns a <i>k</i>-by-1 vector; the function's input
   *           is <b>x</b> stacked on top of <b>v</b>
   * @param px a Gaussian distribution (in the moment parameterization)
   *           over the <i>n</i>-by-1 vector <b>x</b> (which can be 
   *           partitioned into several sub-variables)
   * @throws IllegalArgumentException if the sum dimension of 
   *                                  <code>px</code> and
   *                                  <code>pv</code> does not match the
   *                                  input dimension of <code>f</code> */
  public ExtendedTransformation(NoisyVectorFunction f, Gaussian px) {
    super(f, px);
    ExtendedVectorFunction g = (ExtendedVectorFunction)f;
    // Form the product distribution.
    Gaussian pv = f.getNoiseModel();
    q = new Gaussian(px);
    boolean noise = (pv != null);
    if (noise) q.times(pv, true);
    // Get the Jacobian and break it into pieces corresponding to x and v.
    double[] input = q.getMu(null).getColumnPackedCopy();
    Matrix jacobian = new Matrix(((ExtendedVectorFunction)f).jacobian(input));
    Matrix jX = jacobian.getMatrix(0, k - 1, q.getIndices(xSet));
    Matrix jV = null;
    if (noise) jV = jacobian.getMatrix(0, k - 1, q.getIndices(vSet));
    // Compute the parameters of the linear-Gaussian model.
    Matrix muX = px.getMu(xSet);
    Matrix sigmaXX = px.getSigma(xSet, null);
    Matrix sigmaVV = null;
    Matrix muY = new Matrix(f.evaluate(input), y.dim);
    B = jX;
    a = muY.minus(B.times(muX));
    Flops.count(Flops.add(muY.getRowDimension(), 1) +
		Flops.mult(B.getRowDimension(), 
			   B.getColumnDimension(),
			   muX.getRowDimension()));
    if (noise) {
      sigmaVV = pv.getSigma(new ListSet(pv.getVariables()), null);
      G = jV.times(sigmaVV.times(jV.transpose()));
      Flops.count(Flops.mult(jV.getRowDimension(), 
			     jV.getColumnDimension(),
			     sigmaVV.getColumnDimension()) +
		  Flops.mult(jV.getRowDimension(), 
			     sigmaVV.getColumnDimension(),
			     jV.getRowDimension()));
    } else
      G = new Matrix(y.dim, y.dim);
    // Add the function output to the Gaussian distribution q.
    q.extend(ySet);
    q.setMu(ySet, muY);
    q.setSigma(ySet, null, B.times(sigmaXX.times(B.transpose())).plus(G));
    q.setSigma(ySet, xSet, B.times(sigmaXX));
    if (noise) q.setSigma(ySet, vSet, jV.times(sigmaVV));
    Flops.count(Flops.mult(B.getRowDimension(), 
			   B.getColumnDimension(),
			   sigmaXX.getColumnDimension()) +
		Flops.mult(B.getRowDimension(), 
			   sigmaXX.getColumnDimension(),
			   B.getRowDimension()) +
		Flops.mult(jV.getRowDimension(), 
			   jV.getColumnDimension(), 
			   sigmaVV.getColumnDimension()));
  }

  /**
   * Returns a handle on a factory for creating linearizations using
   * the extended transformation.
   *
   * @return a factory for creating linearizations
   */
  public static LinearizationFactory getFactory() {
    return factory;
  }
}
