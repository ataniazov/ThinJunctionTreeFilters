package javaslam.slam;

import java.util.*;
import Jama.*;
import javaslam.util.*;
import javaslam.prob.*;
import javaslam.filter.*;

/**
 * A Kalman filter for the Simultaneous Localization and Mapping (SLAM)
 * problem.
 *
 * <p>This class records counts of all floating point operations using
 * {@link Flops#count(long)} (except those used in the service of
 * debugging and avoiding numerical errors).</p>
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.7 $ ($Date: 2003/02/07 02:26:18 $) 
 */
public class KalmanSLAMFilter extends LGSLAMFilter {

  /**
   * The underlying Kalman filter.
   */
  protected KalmanFilter kf;

  /**
   * For convenience, a list-set containing only {@link SLAMFilter#x}.
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
  public KalmanSLAMFilter(double[] mu, 
			  double[][] sigma) {
    super(mu.length);
    xSet = new ListSet();
    xSet.add(x);
    Gaussian p = new Gaussian(xSet, mu, sigma, true);
    p.setDoubling(true);
    kf = new KalmanFilter(p);
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
   *             linear coefficient for the landmark's state; if 
   *             <code>id</code> identifies a new landmark, then D must 
   *             be invertible
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
    Matrix DM = new Matrix(D);
    Variable lm = getLandmarkVariable(id);
    if (lm == null) {
      lm = addLandmark(id, DM.getColumnDimension());
      Matrix z0M = new Matrix(z0, z0.length);
      Matrix CM = new Matrix(C);
      Matrix RM = new Matrix(R);
      Matrix zM = new Matrix(z, z.length);
      Matrix muX = kf.getDistribution().getMu(xSet);
      Matrix Dinv = DM.inverse();
      Flops.count(Flops.inv(DM.getColumnDimension(), false));
      ListSet lmSet = new ListSet(lm);
      kf.getDistribution().extend(lmSet);
      kf.getDistribution().
	setMu(lmSet, Dinv.times(zM.minus(z0M).minus(CM.times(muX))));
      Flops.count(Flops.mult(CM.getRowDimension(), 
			     CM.getColumnDimension(), 
			     muX.getColumnDimension()) +
		  2 * Flops.add(zM.getRowDimension(), 
				zM.getColumnDimension()) +
		  Flops.mult(DM.getRowDimension(), 
			     DM.getColumnDimension(), 
			     zM.getColumnDimension()));
      // Compute the landmark state covariance.
      Matrix sigmaXX = kf.getDistribution().getSigma(xSet, null);
      Matrix sigmaLX = CM.times(sigmaXX);
      Matrix tmp = sigmaLX.times(CM.transpose()).plus(RM);
      Matrix sigmaLL = Dinv.times(tmp).times(Dinv.transpose());
      kf.getDistribution().setSigma(lmSet, null, sigmaLL);
      Flops.count(Flops.mult(CM.getRowDimension(), 
			     CM.getColumnDimension(), 
			     sigmaXX.getColumnDimension()) +
		  Flops.mult(sigmaLX.getRowDimension(), 
			     sigmaLX.getColumnDimension(), 
			     CM.getRowDimension()) +
		  Flops.add(RM.getRowDimension(), 
			    RM.getColumnDimension()) +
		  Flops.mult(Dinv.getRowDimension(), 
			     Dinv.getColumnDimension(), 
			     tmp.getColumnDimension()) +
		  Flops.mult(Dinv.getRowDimension(), 
			     tmp.getColumnDimension(), 
			     Dinv.getRowDimension()));
      // Compute the covariances between the landmark state and the
      // robot and other landmark states.
      Matrix tmp2 = Dinv.times(CM).transpose().uminus();
      Matrix sigmaXL = sigmaXX.times(tmp2);
      kf.getDistribution().setSigma(xSet, lmSet, sigmaXL);
      ListSet lmsSet = new ListSet(id2lm.values());
      lmsSet.remove(lm);
      Matrix sigmaUX = kf.getDistribution().getSigma(lmsSet, xSet);
      Matrix sigmaUL = sigmaUX.times(tmp2);
      kf.getDistribution().setSigma(lmsSet, lmSet, sigmaUL);
      // Count flops.
      Flops.count(Flops.mult(Dinv.getRowDimension(),
			     Dinv.getColumnDimension(), 
			     CM.getColumnDimension()) +
		  Flops.mult(sigmaXX.getRowDimension(), 
			     sigmaXX.getColumnDimension(), 
			     tmp2.getColumnDimension()) +
		  Flops.mult(sigmaUX.getRowDimension(), 
			     sigmaUX.getColumnDimension(), 
			     tmp2.getColumnDimension()));
      return;
    }
    ListSet vars = new ListSet(x);
    vars.add(lm);
    Matrix E = new Matrix(z.length, x.dim + lm.dim);
    E.setMatrix(0, z.length - 1, 0, x.dim - 1, new Matrix(C));
    E.setMatrix(0, z.length - 1, x.dim, x.dim + lm.dim - 1, DM);
    kf.measurement(vars, z0, E.getArray(), R, z);
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
    kf.time(xSet, x0, A, Q);
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
  public void odometry(double[] y0, 
		       double[][] B, 
		       double[][] S, 
		       double[] y) {
    kf.measurement(new ListSet(x), y0, B, S, y);
  }

  /**
   * Returns the filtered marginal potential over the robot's state.
   *
   * @return the filtered marginal potential over the robot's state.
   */
  public Gaussian getRobotMarginal() {
    return kf.getMarginal(xSet);
  }

  /**
   * Returns the filtered marginal potential over a landmark's state
   * (in the moment parameterization).
   *
   * @param ids the landmark identifier
   * @return the filtered marginal potential over a landmark's state
   */
  public Gaussian getLandmarkMarginal(int id) {
    return kf.getMarginal(Collections.singleton(getLandmarkVariable(id)));
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
    return kf.getMarginal(tmp);
  }

  /**
   * Returns the underlying Kalman filter used by this SLAM filter.
   *
   * @return the underlying Kalman filter used by this SLAM filter
   */
  public KalmanFilter getKalmanFilter() {
    return kf;
  }
}
