package javaslam.filter;

import java.util.*;
import Jama.*;
import javaslam.prob.*;
import javaslam.util.*;

/**
 * An Information filter.
 *
 * <p>This class records counts of all floating point operations using
 * {@link Flops#count(long)} (except those used in the service of
 * debugging and avoiding numerical errors).</p>
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.2 $ ($Date: 2003/02/07 20:06:18 $)
 */
public class InformationFilter implements Filter {

  /**
   * The canonical potential representing the belief state.
   */
  protected Gaussian p;

  /**
   * Constructor.
   *
   * @param p the initial belief state
   */
  public InformationFilter(Gaussian p) {
    if (p.isMoment())
      p.reparameterize(true);
    this.p = p;
  }

  /**
   * Returns the current filtered belief state.
   *
   * @return the current filtered belief state
   */
  public Gaussian getDistribution() {
    return p;
  }

  /**
   * Marginalizes a set of variables out of the belief state.
   *
   * @param mvars a set of {@link Variable Variable}s to marginalize out
   */
  public void marginalizeOut(Set mvars) {
    p.marginalizeOut(mvars);
  }

  /**
   * Creates a linear-Gaussian measurement potential.  The parameters
   * <code>vars</code>, <code>C</code> and <code>R</code> define the
   * measurement equation as follows:
   *
   * <blockquote>
   * <nobr><b>y</b>=<code>y0</code> + 
   * <code>C</code><b>x</b>(<code>vars</code>)+<b>w</b></nobr>
   * </blockquote>
   *
   * where <b>w</b> is a white-noise variable with covariance
   * <code>R</code>.  Given the actual measurement <code>y</code>,
   * this method computes the potential over <code>vars</code> that,
   * when multiplied into a distribution over <code>vars</code>, has
   * the effect of conditioning on <code>y</code>.
   *
   * @param vars an ordered set of the variables with sum dimension 
   *             <i>n</i> in the belief state that causally influenced
   *             this measurement; any variables in this list that are
   *             not currently in the belief state are added with
   *             uninformative priors.
   * @param y0   a <i>k</i>-vector giving the constant term
   * @param C    a <i>k</i> by <i>n</i> observation matrix that defines the
   *             linear measurement model (and whose columns are ordered 
   *             consistently with the order of <code>vars</code>)
   * @param R    a <i>k</i> by <i>k</i> symmetric positive definite matrix
   *             giving the covariance of the measurement white noise
   * @param y    the measurement <i>k</i>-vector
   * @return a Gaussian in the canonical parameterization which, when
   *         multiplied into a graphical model, has the effect of 
   *         conditioning on the observation <code>y</code>
   * @throws IllegalArgumentException if there are any dimension mismatches
   */
  public static Gaussian getConditionalGaussian(ListSet vars, 
						Matrix y0, Matrix C, 
						Matrix R, Matrix y) {
    // Check dimensions.
    int dim = Variable.dimension(vars);
    if (y0.getColumnDimension() != y.getColumnDimension())
      throw new IllegalArgumentException("The length of y0 does " +
					 "not match the length of y.");
    if (C.getColumnDimension() != dim)
      throw new IllegalArgumentException("The number of columns of C does " +
					 "not match the sum dimension of " +
					 "the variables.");
    if (C.getRowDimension() != y.getRowDimension())
      throw new IllegalArgumentException("The number of rows of C does " +
					 "not match the length of y.");
    if (R.getRowDimension() != R.getColumnDimension())
      throw new IllegalArgumentException("R is not square.");
    if (R.getRowDimension() != y.getRowDimension())
      throw new IllegalArgumentException("The number of rows of R does " +
					 "not match the length of y.");
    // Construct the potential.
    Matrix t = R.solve(C);
    Gaussian p = new Gaussian(vars, t.transpose().times(y.minus(y0)),
				C.transpose().times(t), false);
    // Count flops.
    Flops.count(Flops.solve(R.getRowDimension(), 
			    R.getColumnDimension(), 
			    C.getColumnDimension(), 
			    true) +
		Flops.add(y.getRowDimension(), 
			  y.getColumnDimension()) +
		Flops.mult(t.getColumnDimension(), 
			   t.getRowDimension(), 
			   y.getColumnDimension()) +
		Flops.mult(C.getColumnDimension(), 
			   C.getRowDimension(), 
			   t.getColumnDimension()));
    return p;
  }

  /**
   * Performs a linear-Gaussian measurement update.  The parameters
   * <code>vars</code>, <code>C</code> and <code>R</code> define the
   * measurement equation as follows:
   *
   * <blockquote>
   * <nobr><b>y</b>=<code>y0</code> + 
   * <code>C</code><b>x</b>(<code>vars</code>)+<b>w</b></nobr>
   * </blockquote>
   *
   * where <b>w</b> is a white-noise variable with covariance
   * <code>R</code>.  Given the actual measurement <code>y</code>,
   * this method updates the belief state in <i>O</i>(1) time.
   *
   * @param vars an ordered set of the variables with sum dimension 
   *             <i>n</i> in the belief state that causally influenced
   *             this measurement; any variables in this list that are
   *             not currently in the belief state are added with
   *             uninformative priors.
   * @param y0   a <i>k</i>-vector giving the constant term
   * @param C    a <i>k</i> by <i>n</i> observation matrix that defines the
   *             linear measurement model (and whose columns are ordered 
   *             consistently with the order of <code>vars</code>)
   * @param R    a <i>k</i> by <i>k</i> symmetric positive definite matrix
   *             giving the covariance of the measurement white noise
   * @param y    the measurement <i>k</i>-vector
   * @throws IllegalArgumentException if there are any dimension mismatches */
  public void measurement(ListSet vars, double[] y0, double[][] C, 
			  double[][] R, double[] y) {
    Gaussian q = getConditionalGaussian(vars, 
					new Matrix(y0, y0.length), 
					new Matrix(C), 
					new Matrix(R), 
					new Matrix(y, y.length));
    p.times(q, true);
  }

  /**
   * Constructs a linear-Gaussian potential that introduces a set of
   * new unobserved variables.  The parameters <code>vars</code>,
   * <code>newVars</code>, <code>x0</code>, <code>A</code> and
   * <code>Q</code> define the distribution of the new variables as
   * follows:
   *
   * <blockquote>
   * <nobr><code>newVars</code>=<code>x0</code>
   * <code>A</code><code>vars</code>+<b>v</b></nobr>
   * </blockquote>
   *
   * where <b>v</b> is a white-noise variable with covariance
   * <code>Q</code>.
   *
   * <p>Multiplying the returned potential into a junction tree will not
   * result in message passing; the junction tree is already
   * consistent.  This is because the new variables are barren nodes,
   * i.e., unobserved children of the other variables, and thus do not
   * change their distribution.</p>
   *
   * @param newVars an ordered set of variables with sum dimension 
   *                <i>n</i> that are not present in the belief state
   * @param vars    an ordered set of variables with sum dimension 
   *                <i>m</i> that are currently in the belief state
   * @param x0      an <i>n</i>-vector giving the constant term
   * @param A       an <i>n</i> by <i>m</i> linear coefficient matrix
   * @param Q       an <i>n</i> by <i>n</i> symmetric positive definite 
   *                matrix giving the covariance of the white noise
   * @return a Gaussian in the canonical parameterization which, when
   *         multiplied into a graphical model, has the effect of 
   *         adding <code>newVars</code> into the belief state as a 
   *         directed child of <code>vars</code>
   * @throws IllegalArgumentException if there are any dimension mismatches
   */
  public static Gaussian getUnconditionalGaussian(ListSet newVars, 
						  ListSet vars, 
						  Matrix x0, 
						  Matrix A, 
						  Matrix Q) {
    // Multiply in the potential
    int m = Variable.dimension(vars);
    int n = Variable.dimension(newVars);
    // Check dimensions.
    if (x0.getRowDimension() != n)
      throw new IllegalArgumentException("The length of x0 does " +
					 "not match the sum dimension of " +
					 "the variables");
    if (A.getRowDimension() != n)
      throw new IllegalArgumentException("The number of rows of A does " +
					 "not match the sum dimension of " +
					 "the new variables");
    if (A.getColumnDimension() != m)
      throw new IllegalArgumentException("The number of columns of A does " +
					 "not match the sum dimension of " +
					 "the current variables");
    if (Q.getRowDimension() != Q.getColumnDimension())
      throw new IllegalArgumentException("Q is not square.");
    if (Q.getRowDimension() != n)
      throw new IllegalArgumentException("The number of rows of Q does " +
					 "not match the sum dimension of " +
					 "the new variables");
    // Make the potential over all the variables.
    ListSet allVars = new ListSet();
    allVars.addAll(newVars);
    allVars.addAll(vars);
    /* Since we've stacked our variables like
     *
     *     [ newVars;
     *       vars      ]
     *
     * the parameters of the potential are
     *
     * n = [  Q^{-1} x0;
     *       -A' Q^{-1} x0; ]
     *
     * L = [ Q^{-1},    -A' Q^{-1}; 
     *       -Q^{-1} A,  A' Q^{-1} A ]
     */
    Gaussian p = new Gaussian(allVars, false);
    // Regularize the noise covariance matrix.
    TJTF.regularize(Q, 1e5);
    Matrix Qinv = Q.inverse();
    Matrix t = A.transpose().times(Qinv);
    p.setEta(vars, t.times(x0).uminus());
    p.setEta(newVars, Qinv.times(x0));
    p.setLambda(newVars, null, Qinv);
    p.setLambda(vars, null, t.times(A));
    p.setLambda(vars, newVars, t.uminus());
    // Count flops.
    Flops.count(Flops.inv(Q.getRowDimension(), true) +
		Flops.mult(A.getColumnDimension(), 
			   A.getRowDimension(), 
			   Q.getColumnDimension()) +
		Flops.mult(t.getRowDimension(), 
			   t.getColumnDimension(), 
			   x0.getColumnDimension()) +
		Flops.mult(Qinv.getRowDimension(), 
			   Qinv.getColumnDimension(), 
			   x0.getColumnDimension()) +
		Flops.mult(t.getRowDimension(), 
			   t.getColumnDimension(), 
			   A.getColumnDimension()));
    return p;
  }

  /**
   * Performs a linear-Gaussian time update.  The parameters
   * <code>vars</code>, <code>A</code> and <code>Q</code> define the
   * state evolution equation as follows:
   *
   * <blockquote>
   * <nobr><b>x</b><sub>t + 1</sub>(<code>vars</code>)=<code>x0</code>
   * <code>A</code><b>x</b><sub>t</sub>(<code>vars</code>)+<b>v</b></nobr>
   * </blockquote>
   *
   * where <b>v</b> is a white-noise variable with covariance
   * <code>Q</code>.  All variables not in <code>vars</code> are
   * assumed stationary.
   *
   * @param vars an ordered set of the variables with sum dimension <i>n</i>
   *             in the belief state that evolve over time
   * @param x0   an <i>n</i>-vector giving the constant term
   * @param A    an <i>n</i> by <i>n</i> evolution matrix that defines the
   *             linear evolution model (and whose blocks are ordered 
   *             consistently with the order of <code>vars</code>)
   * @param Q    an <i>n</i> by <i>n</i> symmetric positive definite matrix
   *             giving the covariance of the evolution white noise (and 
   *             whose blocks are ordered consistently with the order
   *             of <code>vars</code>)
   * @throws IllegalArgumentException if there are any dimension mismatches 
   *                                  or <code>vars</code> contains variables
   *                                  that are not in the current belief state
   */
  public void time(ListSet vars, double[] x0, 
		   double[][] A, double[][] Q) {
    // Create anonymous variables that represent vars at the next time step.
    ListSet nextVars = new ListSet();
    for (Iterator i = vars.iterator(); i.hasNext();) {
      Variable v = (Variable)i.next();
      nextVars.add(new Variable(v.label + "\'", v.dim));
      if (!p.getVariables().contains(v))
	throw new IllegalArgumentException("Variable " + v + 
					   " is not in the belief state.");
    }
    Gaussian q = getUnconditionalGaussian(nextVars, vars, 
					  new Matrix(x0, x0.length), 
					  new Matrix(A), new Matrix(Q));
    p.times(q, true);
    // Marginalize out the variables at the previous time step.
    p.marginalizeOut(vars);
    // Now substitute them for their anonymous variables in the model.
    for (Iterator i = vars.iterator(), j = nextVars.iterator();
	 i.hasNext();) {
      Variable v = (Variable)i.next();
      Variable u = (Variable)j.next();
      p.rename(u, v);
    }
  }

  /**
   * Extracts the filtered marginal distribution.
   *
   * @param vars  the set of {@link Variable Variable}s whose filtered
   *              marginal is to be computed
   * @return a filtered marginal potential over <code>vars</code>
   */
  public Gaussian getMarginal(Set vars) {
    return p.marginalize(vars, false);   
  }

  /**
   * Extracts a set of unary marginals.
   *
   * @param  vars a collection of {@link Variable Variable}s, or 
   *         <code>null</code> to indicate all variables in the belief state
   * @return a map whose keys are the (distinct) elements of
   *         <code>vars</code> and whose values are the corresponding 
   *         marginals (in the moment parameterization)
   */
  public Map getMarginals(Collection vars) {
    if (vars == null) vars = getVariables();
    HashMap marginals = new HashMap(vars.size());
    Gaussian q = p.reparameterize(false);
    for (Iterator i = vars.iterator(); i.hasNext();) {
      Variable v = (Variable)i.next();
      marginals.put(v, q.marginalize(Collections.singleton(v), false));
    }
    return marginals;
  }

  /**
   * Gets an unmodifiable set of the {@link Variable Variable}s in the
   * filtered belief state.
   */
  public Set getVariables() {
    return p.getVariables();
  }
}
