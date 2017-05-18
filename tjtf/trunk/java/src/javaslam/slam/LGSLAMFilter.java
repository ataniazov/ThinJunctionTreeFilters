package javaslam.slam;

import java.util.*;
import javaslam.prob.*;

/**
 * A linear-Gaussian filter for the Simultaneous Localization and
 * Mapping (SLAM) problem.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.4 $ ($Date: 2003/02/07 02:26:18 $) */
public abstract class LGSLAMFilter extends AbstractSLAMFilter {

  /**
   * Constructor.
   *
   * @param dim the dimension of the robot state variable
   */
  public LGSLAMFilter(int dim) {
    super(dim);
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
  public abstract void measurement(int id, 
				   double[] z0, 
				   double[][] C, 
				   double[][] D, 
				   double[][] R, 
				   double[] z);

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
  public abstract void motion(double[] x0, double[][] A, double[][] Q);

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
  public abstract void odometry(double[] y0, 
				double[][] B, 
				double[][] S, 
				double[] y);

  /**
   * Returns the filtered marginal potential over the robot's state
   * (in the moment parameterization).
   *
   * @return the filtered marginal potential over the robot's state.
   */
  public abstract Gaussian getRobotMarginal();

  /**
   * Returns the filtered marginal potential over a landmark's state
   * (in the moment parameterization).
   *
   * @param ids the landmark identifier
   * @return the filtered marginal potential over a landmark's state
   */
  public abstract Gaussian getLandmarkMarginal(int id);

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
    Gaussian[] p = new Gaussian[ids.length];
    for (int i = 0; i < ids.length; i++)
      p[i] = getLandmarkMarginal(ids[i]);
    return p;
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
  public abstract Gaussian getRobotLandmarkMarginal(int id);

  /**
   * Given a set of landmarks { <i>l</i><sub>1</sub>, ...,
   * <i>l</i><sub>k</sub>}, this method computes a set of potentials
   * over (<i>x</i>, <i>l</i><sub>1</sub>), ..., (<i>x</i>,
   * <i>l</i><sub>k</sub>).  Each potential is either the filtered
   * marginal over the robot and landmark states.
   *
   * @param ids an array of landmark identifiers (or
   *            <code>null</code> to indicate all landmarks)
   * @return an array of (perhaps approximate) marginal potentials
   *         over the corresponding landmark and robot state 
   *         variables (in the moment parameterization)
   */
  public Gaussian[] getRobotLandmarkMarginals(int[] ids) {
    if (ids == null) ids = getLandmarkIds();
    Gaussian[] p = new Gaussian[ids.length];
    for (int i = 0; i < ids.length; i++)
      p[i] = getRobotLandmarkMarginal(ids[i]);
    return p;
  }

  /**
   * Returns the filtered estimate of the robot's state.
   *
   * @return the filtered estimate of the robot's state
   */
  public double[] getRobotEstimate() {
    return getRobotMarginal().getMu(null).getColumnPackedCopy();
  }

  /**
   * Returns the filtered estimate of a landmark's state.
   *
   * @param id the identifier of the landmark
   * @return the filtered estimate of the landmark's state
   */
  public double[] getLandmarkEstimate(int id) {
    return getLandmarkMarginal(id).getMu(null).getColumnPackedCopy();
  }
}
