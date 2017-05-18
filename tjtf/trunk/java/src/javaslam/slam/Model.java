package javaslam.slam;

import javaslam.filter.*;

/**
 * Represents a <font size="-1">SLAM</font> model.
 * 
 * @author Mark A. Paskin
 * @version $Revision: 1.1 $ ($Date: 2003/02/07 02:26:18 $) 
 */
public interface Model {

  /**
   * Returns the dimension of the robot's state.
   *
   * @return the dimension of the robot's state
   */
  public int getRobotDim();

  /**
   * Returns the dimension of each landmark's state.
   *
   * @return the dimension of each landmark's state.
   */
  public int getLandmarkDim();

  /**
   * Returns the motion model <nobr><b>x</b><sub><i>t</i> + 1</sub> =
   * <i>f</i>(<b>x</b><sub><i>t</i></sub>, <b>w</b>)</nobr>.  In other
   * words, this returns the function that computes the state
   * <b>x</b><sub><i>t</i></sub> the robot would be in given its state is
   * <b>x</b><sub><i>t</i></sub>, it applied the control
   * <code>c</code>, and the motion noise was <b>w</b>.
   *
   * @return the motion model 
   */
  public NoisyVectorFunction getMotionModel(double[] c);

  /**
   * Returns the odometry model <nobr><b>y</b><sub><i>t</i></sub> =
   * <i>g</i>(<b>x</b><sub><i>t</i></sub>, <b>v</b>)</nobr>.  In other
   * words, this returns the function that computes the odometry
   * measurement <b>y</b><sub><i>t</i></sub> the robot would receive
   * given its state is <b>x</b><sub><i>t</i></sub> and the
   * observation noise is <b>v</b>.
   *
   * @return the odometry model 
   */
  public NoisyVectorFunction getOdometryModel();

  /**
   * Returns the landmark observation model
   * <nobr><b>z</b><sub><i>t</i></sub><sup><i>i</i></sup> =
   * <i>h</i>(<b>x</b><sub><i>t</i></sub>,
   * <b>l</b><sub><i>j</i></sub>, <b>u</b>)</nobr>.  In other words,
   * this returns the function that computes the observation
   * <b>z</b><sub><i>t</i></sub><sup><i>i</i></sup> the robot would
   * receive of landmark <i>j</i> given its state is
   * <b>x</b><sub><i>t</i></sub>, the state of the landmark is
   * <b>l</b><sub><i>j</i></sub>, and the observation noise was
   * <b>u</b>.
   *
   * @return the landmark observation model 
   */
  public NoisyVectorFunction getMeasurementModel();

  /**
   * Returns the inverse of the landmark observation model (at
   * <b>z</b><sub><i>t</i></sub><sup><i>i</i></sup> = <code>z</code>)
   * with respect to the landmark state:
   * <nobr><b>l</b><sub><i>j</i></sub> = <i>h</i><sup>-1</sup><sub>
   * <code>z</code></sub>(<b>x</b><sub><i>t</i></sub>,
   * <b>u</b>)</nobr>.  In other words, this returns the function
   * that, given the state of the robot <b>x</b><sub><i>t</i></sub>
   * and the observation noise <b>u</b>, computes the state of a
   * landmark that yielded the observation <code>z</code>.
   *
   * @param z the landmark measurement */
  public NoisyVectorFunction getInverseMeasurementModel(double[] z);
}
