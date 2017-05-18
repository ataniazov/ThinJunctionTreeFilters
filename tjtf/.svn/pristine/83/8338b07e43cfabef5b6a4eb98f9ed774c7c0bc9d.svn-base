package javaslam.slam;

import java.util.*;
import Jama.*;
import javaslam.util.*;
import javaslam.prob.*;
import javaslam.filter.*;
import javaslam.tjt.*;
import javaslam.tjt.graph.*;

/**
 * An information filter for the Simultaneous Localization and
 * Mapping (SLAM) problem.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.3 $ ($Date: 2003/02/07 20:06:31 $) 
 */
public class InformationSLAMFilter extends LGSLAMFilter {

  /**
   * The thin junction tree used to represent the filtered belief
   * state.
   */
  protected InformationFilter pf;

  /**
   * For convenience, a list-set containing only {@link #x}.
   */
  private final ListSet xSet;

  /**
   * Constructor.
   *
   * @param mu      a <i>n</i>-vector containing the starting state 
   *                of the robot
   * @param sigma   an <i>n</i>-by-<i>n</i> positive semidefinite 
   *                covariance matrix giving the uncertainty in the 
   *                robot's initial state
   */
  public InformationSLAMFilter(double[] mu, 
			       double[][] sigma) {
    super(mu.length);
    xSet = new ListSet();
    xSet.add(x);
    Gaussian p = new Gaussian(xSet, mu, sigma, true);
    p.setDoubling(true);
    pf = new InformationFilter(p);
  }

  /**
   * Performs a linear-Gaussian landmark measurement update.  The
   * parameters <code>z0</code>, <code>C</code>, <code>D</code>, and
   * <code>R</code> define the measurement equation as follows:
   *
   * <blockquote>
   * <nobr><b>z</b>=<code>z0</code> + 
   * <code>C</code><b>x</b> + <code>D</code><code>lm</code> + <b>w</b></nobr>
   * </blockquote>
   *
   * where <b>x</b> is the robot state variable, <code>lm</code> is a
   * landmark variable, and <b>w</b> is a white-noise variable with
   * covariance <code>R</code>.
   *
   * @param id   the identifier of the landmark that generated this 
   *             observation; if it is not currently in the belief state, 
   *             then it is first added with an uninformative prior.
   * @param z0   a <i>k</i>-vector giving the constant term
   * @param C    a <i>k</i> by <i>n</i> matrix that defines the
   *             linear coefficient for the robot's state
   * @param D    a <i>k</i> by <i>m</i> matrix that defines the
   *             linear coefficient for the landmark's state
   * @param R    a <i>k</i> by <i>k</i> symmetric positive definite matrix
   *             giving the covariance of the measurement white noise
   * @param z    the measurement <i>k</i>-vector
   * @throws IllegalArgumentException if there are any dimension mismatches 
   */
  public void measurement(int id, 
			  double[] z0, 
			  double[][] C, 
			  double[][] D, 
			  double[][] R, 
			  double[] z) {
    // Construct the measurement potential over the robot's and
    // landmark's state variables.
    Matrix DM = new Matrix(D);
    Variable lm = getLandmarkVariable(id);
    if (lm == null) lm = addLandmark(id, DM.getColumnDimension());
    ListSet vars = new ListSet();
    vars.add(x);
    vars.add(lm);
    Matrix CD = new Matrix(z0.length, x.dim + lm.dim);
    CD.setMatrix(0, z0.length - 1, 0, x.dim - 1, new Matrix(C));
    CD.setMatrix(0, z0.length - 1, x.dim, x.dim + lm.dim - 1, DM);
    pf.measurement(vars, z0, CD.getArray(), R, z);
  }

  /**
   * Performs a linear-Gaussian motion update.  The parameters
   * <code>A</code> and <code>Q</code> define the state evolution
   * equation as follows:
   *
   * <blockquote>
   * <nobr><b>x</b><sub>t + 1</sub> = <code>x0</code> +
   * <code>A</code><b>x</b><sub>t</sub>+<b>v</b></nobr>
   * </blockquote>
   *
   * where <b>v</b> is a white-noise variable with covariance
   * <code>Q</code>.
   *
   * @param x0   an <i>n</i>-vector giving the constant term of the 
   *             state evolution equation
   * @param A    an <i>n</i> by <i>n</i> evolution matrix that defines the
   *             linear evolution model
   * @param Q    an <i>n</i> by <i>n</i> symmetric positive definite matrix
   *             giving the covariance of the evolution white noise
   * @throws IllegalArgumentException if there are any dimension mismatches 
   */
  public void motion(double[] x0, double[][] A, double[][] Q) {
    pf.time(xSet, x0, A, Q);
  }

  /**
   * Performs a linear-Gaussian odometry update.  The parameters
   * <code>y0</code>, <code>B</code>, and <code>S</code> define the
   * odometry measurement equation as follows:
   *
   * <blockquote>
   * <nobr><b>y</b><sub>t + 1</sub>=<code>y0</code> +
   * <code>B</code><b>x</b><sub>t + 1</sub>+<b>u</b></nobr>
   * </blockquote>
   *
   * where <b>u</b> is a white-noise variable with covariance
   * <code>S</code>.
   *
   * @param y0   an <i>h</i>-vector giving the constant term of the 
   *             odometry equation
   * @param B    an <i>h</i> by <i>n</i> evolution matrix that is
   *             linear term of the odometry equation
   * @param S    an <i>h</i> by <i>h</i> symmetric positive definite matrix
   *             giving the covariance of the odometry white noise
   * @param y    an <i>h</i>-vector giving the odometry measurement
   * @throws IllegalArgumentException if there are any dimension mismatches 
   */
  public void odometry(double[] y0, double[][] B, double[][] S, double[] y) {
    pf.measurement(xSet, y0, B, S, y);
  }

  /**
   * Returns the filtered marginal potential (in the moment
   * parameterization) over the robot's state.
   *
   * @return the filtered marginal potential over the robot's state.
   */
  public Gaussian getRobotMarginal() {
    return pf.getMarginal(xSet).reparameterize(true);
  }

  /**
   * Returns the filtered marginal potential over a landmark's state
   * (in the moment parameterization).
   *
   * @param ids the landmark identifier
   * @return the filtered marginal potential over a landmark's state
   */
  public Gaussian getLandmarkMarginal(int id) {
    return pf.getMarginal(Collections.singleton(getLandmarkVariable(id))).
      reparameterize(true);
  }

  /**
   * This method computes the joint filtered distribution over the
   * states of the robot and a landmark.
   *
   * @param id the landmark's identifier
   * @return the joint filtered distribution over the
   *         states of the robot and a landmark (in the moment 
   *         parameterization)
   */
  public Gaussian getRobotLandmarkMarginal(int id) {
    ListSet tmp = new ListSet(getRobotVariable());
    tmp.add(getLandmarkVariable(id));
    return pf.getMarginal(tmp).reparameterize(true);
  }

  /**
   * Given a collection of landmark state variables, this method
   * returns the set of filtered marginal potentials over each
   * individual landmark's state.
   *
   * @param ids an array of landmark identifiers (or
   *            <code>null</code> to indicate all landmarks)
   * @return a corresponding array of marginals in the moment
   *         parameterization
   */
  public Gaussian[] getLandmarkMarginals(int[] ids) {
    if (ids == null) ids = getLandmarkIds();
    Gaussian[] m = new Gaussian[ids.length];
    Gaussian q = pf.getDistribution().reparameterize(false);
    for (int i = 0; i < ids.length; i++) {
      Variable landmark = getLandmarkVariable(ids[i]);
      m[i] = q.marginalize(Collections.singleton(landmark), false);
    }
    return m;
  }

  /**
   * Given a set of landmarks { <i>l</i><sub>1</sub>, ...,
   * <i>l</i><sub>k</sub>}, this method computes a set of potentials
   * over (<i>x</i>, <i>l</i><sub>1</sub>), ..., (<i>x</i>,
   * <i>l</i><sub>k</sub>).  Each potential is either the filtered
   * marginal over the robot and landmark states, or a
   * product-of-marginals approximation thereof.
   *
   * @param ids an array of landmark identifiers (or
   *            <code>null</code> to indicate all landmarks)
   * @return an array of (perhaps approximate) marginal potentials
   *         over the corresponding landmark and robot state 
   *         variables (in the moment parameterization)
   */
  public Gaussian[] getRobotLandmarkMarginals(int[] ids) {
    if (ids == null) ids = getLandmarkIds();
    Gaussian[] m = new Gaussian[ids.length];
    Gaussian q = pf.getDistribution().reparameterize(false);
    ListSet tmp = new ListSet(getRobotVariable());
    for (int i = 0; i < ids.length; i++) {
      Variable landmark = getLandmarkVariable(ids[i]);
      tmp.add(landmark);
      m[i] = q.marginalize(tmp, false);
      tmp.remove(landmark);
    }
    return m;
  }
}
