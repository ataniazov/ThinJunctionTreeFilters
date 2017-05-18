package javaslam.slam;

import javaslam.filter.*;

/**
 * <p>An interface implemented by classes that implement nonlinear
 * filters for the Simultaneous Localization and Mapping (SLAM)
 * problem.</p>
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.1 $ ($Date: 2003/02/07 02:26:18 $) 
 */
public interface NonlinearSLAMFilter extends SLAMFilter {

  /**
   * Performs a motion update.
   *
   * @param f the stochastic motion model
   */
  public abstract void motion(NoisyVectorFunction f);

  /**
   * Performs an odometry update.
   *
   * @param g the noisy odometry measurement model
   * @param y the odometry measurement
   */
  public abstract void odometry(NoisyVectorFunction g, double[] y);

  /**
   * Performs a landmark measurement update assuming a known data
   * association.
   *
   * @param h    the stochastic landmark measurement model
   * @param hInv the inverse of h with respect to the landmark state
   * @param id   the identifier of the landmark that was observed
   * @param z    the landmark measurement
   */
  public abstract void measurement(NoisyVectorFunction h,  
				   NoisyVectorFunction hInv,
				   int id, double[] z);
}
