package javaslam.filter;

import java.util.*;
import Jama.*;
import javaslam.prob.*;
import javaslam.util.*;
import javaslam.tjt.*;

/**
 * A junction tree filter.
 *
 * <p>This class records counts of all floating point operations using
 * {@link Flops#count(long)} (except those used in the service of
 * debugging and avoiding numerical errors).</p>
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.6 $ ($Date: 2003/01/07 01:16:11 $)
 */
public class JTFilter implements Filter {

  /**
   * The junction tree used to represent the belief state.
   */
  protected JunctionTree jt;

  /**
   * Default constructor.
   */
  public JTFilter(JunctionTree jt) {
    this.jt = jt;
  }

  /**
   * Marginalizes a set of variables out of the belief state.
   *
   * @param mvars a set of {@link Variable Variable}s to marginalize out
   */
  public void marginalizeOut(Set mvars) {
    for (Iterator i = mvars.iterator(); i.hasNext();)
      jt.marginalizeOut((Variable)i.next());
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
   * this method updates the belief state.
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
   * @throws IllegalArgumentException if there are any dimension mismatches
   */
  public void measurement(ListSet vars, double[] y0, double[][] C, 
			  double[][] R, double[] y) {
    Gaussian p = InformationFilter.
      getConditionalGaussian(vars, 
			     new Matrix(y0, y0.length), 
			     new Matrix(C), 
			     new Matrix(R), 
			     new Matrix(y, y.length));
    if (TJTF.blather) TJTF.say("Multiplying in " + p);
    jt.times(p);
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
  public void time(ListSet vars, double[] x0, double[][] A, double[][] Q) {
    // Create anonymous variables that represent vars at the next time step.
    ListSet nextVars = new ListSet();
    for (Iterator i = vars.iterator(); i.hasNext();) {
      Variable v = (Variable)i.next();
      nextVars.add(new Variable(v.label + "\'", v.dim));
      if (!jt.contains(v))
	throw new IllegalArgumentException("Variable " + v + 
					   " is not in the belief state.");
    }
    Gaussian p = InformationFilter.
      getUnconditionalGaussian(nextVars, vars, 
			       new Matrix(x0, x0.length), 
			       new Matrix(A), new Matrix(Q));
    jt.times(p);
    // Marginalize out the variables at the previous time step.
    for (Iterator i = vars.iterator(); i.hasNext();) {
      Variable v = (Variable)i.next();
      jt.marginalizeOut(v);
    }
    // Now substitute them for their anonymous variables in the model.
    for (Iterator i = vars.iterator(), j = nextVars.iterator();
	 i.hasNext();) {
      Variable v = (Variable)i.next();
      Variable u = (Variable)j.next();
      jt.rename(u, v);
    }
  }

  /**
   * Extracts the filtered marginal distribution.  The correct
   * filtered marginal is returned (up to approximations made by the
   * underlying junction tree); also, a product-of-marginals
   * approximation may be returned if it can be computed more
   * efficiently than the full marginal.
   *
   * @param vars  the set of {@link Variable Variable}s whose filtered
   *              marginal is to be computed
   * @return a filtered marginal potential over <code>vars</code>
   * @see JunctionTree#getMarginal(Set,boolean)
   */
  public Gaussian getMarginal(Set vars) {
    Gaussian p = jt.getMarginal(vars, false);
    if (p == null) {
      p = new Gaussian(false);
      for (Iterator i = vars.iterator(); i.hasNext();) {
	Gaussian q = jt.getMarginal(Collections.singleton(i.next()), 
				     false);
	p.times(q, true);
      }
    }
    return p;
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
    return jt.getMarginals(vars);
  }

  /**
   * Gets an unmodifiable set of the {@link Variable Variable}s in the
   * filtered belief state.  The iteration order of this set is the order
   * in which the variables were added to the belief state.
   */
  public Set getVariables() {
    return jt.getVariables();
  }

  /**
   * Gets the junction tree representation of the belief state.
   */
  public JunctionTree getJunctionTree() {
    return jt;
  }

  /**
   * Returns a string representation of this filter.
   */
  public String toString() {
    return "Filter with belief state represented as a " + jt.toString();
  }
}
