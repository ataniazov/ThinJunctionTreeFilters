package javaslam.slam;

import javaslam.util.*;
import javaslam.prob.*;
import javaslam.filter.*;

/**
 * A nonlinear filter for the Simultaneous Localization and Mapping
 * <font size="-1">SLAM</font> problem that uses a linear-Gaussian
 * filter with a technique for linearizing nonlinear motion and
 * measurement models.  For example, by coupling a {@link
 * KalmanSLAMFilter} with {@link ExtendedTransformation}
 * linearizations, one obtains the traditional Extended Kalman filter
 * technique for <font size="-1">SLAM</font>.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.1 $ ($Date: 2003/02/07 02:26:18 $) 
 */
public class LinearizedSLAMFilter implements NonlinearSLAMFilter {

  /**
   * The underlying linear-Gaussian SLAM filter.
   */
  protected LGSLAMFilter filter;

  /**
   * The factory for creating linearizations.
   */
  protected LinearizationFactory linearization;

  /**
   * Constructor.
   *
   * @param filter the underlying linear-Gaussian SLAM filter used for
   *               filtering
   * @param linearization a linearization factory used to 
   *                      create linearizations of nonlinear models
   * @see ExtendedTransformation#getFactory()
   * @see UnscentedTransformation#getFactory()
   */
  public LinearizedSLAMFilter(LGSLAMFilter filter,
			      LinearizationFactory linearization) {
    this.filter = filter;
    this.linearization = linearization;
  }

  /**
   * Constructor.  The unscented transformation is used for linearization.
   *
   * @param filter the underlying linear-Gaussian SLAM filter used for
   *               filtering
   * @see javaslam.filter.UnscentedTransformation
   */
  public LinearizedSLAMFilter(LGSLAMFilter filter) {
    this.filter = filter;
    this.linearization = UnscentedTransformation.getFactory();
  }

  /**
   * Performs a motion update.
   *
   * @param f the stochastic motion model
   */
  public void motion(NoisyVectorFunction f) {
    Gaussian p = filter.getRobotMarginal();
    Linearization l = linearization.linearize(f, p);
    filter.motion(l.getConstantTerm().getColumnPackedCopy(),
		  l.getCoefficient().getArray(),
		  l.getNoiseCovariance().getArray());
  }

  /**
   * Performs an odometry update.
   *
   * @param g the noisy odometry measurement model
   * @param y the odometry measurement
   */
  public void odometry(NoisyVectorFunction g, double[] y) {
    Gaussian p = filter.getRobotMarginal();
    Linearization l = linearization.linearize(g, p);
    filter.odometry(l.getConstantTerm().getColumnPackedCopy(),
		    l.getCoefficient().getArray(),
		    l.getNoiseCovariance().getArray(), 
		    y);
  }

  /**
   * Performs a landmark measurement update assuming a known data
   * association.
   *
   * @param h    the stochastic landmark measurement model
   * @param hInv the inverse of h with respect to the landmark state
   * @param id   the identifier of the landmark that was observed
   * @param z    the landmark measurement
   */
  public void measurement(NoisyVectorFunction h, 
			  NoisyVectorFunction hInv, 
			  int id, double[] z) {
    double[]   z0 = null;
    double[][] A = null;
    double[][] B = null;
    double[][] S = null;
    // Get a distribution over the robot and landmark states.
    if (contains(id)) {
      int[] ids = new int[1]; ids[0] = id;
      Gaussian[] p = filter.getRobotLandmarkMarginals(ids);
      Linearization l = linearization.linearize(h, p[0]);
      z0 = l.getConstantTerm().getColumnPackedCopy();
      A = l.getCoefficient(getRobotVariable()).getArray();
      B = l.getCoefficient(getLandmarkVariable(id)).getArray();
      S = l.getNoiseCovariance().getArray();
    } else {
      // If the landmark has not yet been observed, we obtain an
      // approximate Gaussian distribution over its state by
      // linearizing the inverse of the measurement function.
      Gaussian p = filter.getRobotMarginal();
      Linearization il = linearization.linearize(hInv, p);
      ListSet iSet = new ListSet(getRobotVariable());
      iSet.add(il.getOutputVariable());
      Gaussian q = il.getDistribution().marginalize(iSet, true);
      Linearization l = linearization.linearize(h, q);
      z0 = l.getConstantTerm().getColumnPackedCopy();
      A = l.getCoefficient(getRobotVariable()).getArray();
      B = l.getCoefficient(il.getOutputVariable()).getArray();
      S = l.getNoiseCovariance().getArray();
    }
    filter.measurement(id, z0, A, B, S, z);
  }

  /**
   * Returns the underlying linear-Gaussian filter, which can be used
   * to access the filtered belief state.
   *
   * @return the underlying linear-Gaussian filter
   */
  public LGSLAMFilter getLGFilter() {
    return filter;
  }

  /**
   * Gets the state variable associated with the robot.
   */
  public Variable getRobotVariable() {
    return filter.getRobotVariable();
  }

  /**
   * Gets the state variable associated with the landmark with the
   * supplied identifier.
   */
  public Variable getLandmarkVariable(int id) {
    return filter.getLandmarkVariable(id);
  }
  
  /**
   * Gets the identifier associated with the supplied landmark state
   * variable.
   */
  public int getLandmarkId(Variable lm) {
    return filter.getLandmarkId(lm);
  }

  /**
   * Returns the number of landmarks known to this filter.
   */
  public int getNumLandmarks() {
    return filter.getNumLandmarks();
  }

  /**
   * Adds a new landmark to the map and returns its state variable.
   */
  public Variable addLandmark(int id, int dim) {
    return filter.addLandmark(id, dim);
  }

  /**
   * Returns the IDs of landmarks known to this filter.
   */
  public int[] getLandmarkIds() {
    return filter.getLandmarkIds();
  }

  /**
   * Return <code>true</code> if this filter contains the landmark
   * with the supplied identifier.
   */
  public boolean contains(int id) {
    return filter.contains(id);
  }


  /**
   * Returns the filtered estimate of the robot's state.
   *
   * @return the filtered estimate of the robot's state
   */
  public double[] getRobotEstimate() {
    return filter.getRobotEstimate();
  }

  /**
   * Returns the filtered estimate of a landmark's state.
   *
   * @param id the identifier of the landmark
   * @return the filtered estimate of the landmark's state
   */
  public double[] getLandmarkEstimate(int id) {
    return filter.getLandmarkEstimate(id);
  }
}
