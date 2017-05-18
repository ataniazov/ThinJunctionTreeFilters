package javaslam.slam;

import javaslam.util.*;
import javaslam.prob.*;
import javaslam.filter.*;
import Jama.*;

/**
 * A <font size="-1">SLAM</font> model for a robot that can rotate
 * about its vertical axis and translate along its current heading and
 * which receives differential odometry measurements and range-bearing
 * landmark measurements.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.1 $ ($Date: 2003/02/07 02:26:18 $) */
public class ExampleModel implements Model {

  /**
   * The dimension of the robot state.
   */
  protected final static int ROBOT_DIM = 5;
  
  /**
   * The dimension of each landmark's state.
   */
  protected final static int LANDMARK_DIM = 2;

  /**
   * Returns the dimension of the robot's state.
   *
   * @return the dimension of the robot's state
   */
  public int getRobotDim() {
    return ROBOT_DIM;
  }

  /**
   * Returns the dimension of each landmark's state.
   *
   * @return the dimension of each landmark's state.
   */
  public int getLandmarkDim() {
    return LANDMARK_DIM;
  }

  /**
   * A distribution over the noise variables of the motion model.
   */
  protected Gaussian motionNoiseModel;

  /**
   * Constructor.
   */
  public ExampleModel() {
    motionNoiseModel = getMotionNoiseModel();
    odometryNoiseModel = getOdometryNoiseModel();
    landmarkNoiseModel = getLandmarkNoiseModel();
    odometryModel = new OdometryModel();
    landmarkMeasurementModel = new LandmarkMeasurementModel();
  }

  /**
   * Creates a motion noise model.
   *
   * @param tVelRelNoise the variance of the relative noise in the
   *                     translational velocity control signal
   * @param tVelAbsNoise the variance of the absolute noise in the
   *                     translational velocity control signal 
   *                     specified in (meters/second)<sup>2</sup>
   * @param rVelRelNoise the variance of the relative noise in the
   *                     rotational velocity control signal
   * @param rVelAbsNoise the variance of the absolute noise in the
   *                     rotational velocity control signal 
   *                     specified in (radians/second)<sup>2</sup>
   * @return a Gaussian distribution over a vector variable whose
   *         components correspond to the relative and absolute errors 
   *         in the translational and rotational control signals
   */
  protected static Gaussian getMotionNoiseModel(double tVelRelNoise,
						double tVelAbsNoise,
						double rVelRelNoise,
						double rVelAbsNoise) {
    Variable tVelRelNoiseVar = 
      new Variable("Relative noise in translational velocity control", 1);
    Variable tVelAbsNoiseVar = 
      new Variable("Absolute noise in translational velocity control", 1);
    Variable rVelRelNoiseVar = 
      new Variable("Relative noise in rotational velocity control", 1);
    Variable rVelAbsNoiseVar = 
      new Variable("Absolute noise in rotational velocity control", 1);
    ListSet noiseVars = new ListSet(tVelRelNoiseVar);
    noiseVars.add(tVelAbsNoiseVar);
    noiseVars.add(rVelRelNoiseVar);
    noiseVars.add(rVelAbsNoiseVar);
    Matrix mu = new Matrix(4, 1);
    mu.set(0, 0, 1.0d);
    mu.set(1, 0, 0.0d);
    mu.set(2, 0, 1.0d);
    mu.set(3, 0, 0.0d);
    Matrix sigma = new Matrix(4, 4);
    sigma.set(0, 0, tVelRelNoise);
    sigma.set(1, 1, tVelAbsNoise);
    sigma.set(2, 2, rVelRelNoise);
    sigma.set(3, 3, rVelAbsNoise);
    return new Gaussian(noiseVars, mu, sigma, true);
  }

  /**
   * Creates a motion noise model using parameters specified by system
   * properties or default values (when the properties are not defined).
   *
   * @return a Gaussian distribution over a vector variable whose
   *         components correspond to the relative and absolute errors 
   *         in the translational and rotational control signals
   */
  protected static Gaussian getMotionNoiseModel() {
    double tVelRelNoise = 
	Double.parseDouble(System.getProperty("javaslam.CTRL_TVELRELNOISE",
					      Double.toString(0.03 * 0.03)));
    double tVelAbsNoise = 
	Double.parseDouble(System.getProperty("javaslam.CTRL_TVELABSNOISE",
					      Double.toString(0.02 * 0.02)));
    double rVelRelNoise = 
	Double.parseDouble(System.getProperty("javaslam.CTRL_RVELRELNOISE",
					      Double.toString(0.05 * 0.05)));
    double rVelAbsNoise = 
	Double.parseDouble(System.
			   getProperty("javaslam.CTRL_RVELABSNOISE",
				       Double.
				       toString(Math.pow(Math.PI / 180.0d, 2.0))));
    return getMotionNoiseModel(tVelRelNoise,
			       tVelAbsNoise,
			       rVelRelNoise,
			       rVelAbsNoise);
  }

  /**
   * The motion model of the robot.  The control specifies a desired
   * translational and rotational velocity.
   */
  public class MotionModel 
    implements NoisyVectorFunction, ExtendedVectorFunction {

    /**
     * The current translational velocity control (in meters per second).
     */
    protected double tVelCtrl;

    /**
     * The current rotational velocity control (in radians per second).
     */
    protected double rVelCtrl;

    /**
     * The input dimension of the motion model.
     */
    protected final static int INPUT_DIM = 9;

    /**
     * The output dimension of the motion model.
     */
    protected final static int OUTPUT_DIM = 5;
  
    /* These are indexes into the input (and output) vectors that pull
     * out specific elements.  
     */
    protected final static int XPOS = 0;
    protected final static int YPOS = 1;
    protected final static int HEADING = 2;
    protected final static int TVEL = 3;
    protected final static int RVEL = 4;
    protected final static int TVEL_REL_NOISE = 5;
    protected final static int TVEL_ABS_NOISE = 6;
    protected final static int RVEL_REL_NOISE = 7;
    protected final static int RVEL_ABS_NOISE = 8;

    /**
     * Constructor.
     *
     * @param c the control signal; <code>c[0]</code> is the
     *          translational velocity (in meters/second) and
     *          <code>c[1]</code> is the rotational velocity (in
     *          radians/second) 
     */
    public MotionModel(double[] c) {
      tVelCtrl = c[0];
      rVelCtrl = c[1];
    }

    /**
     * Returns the input dimension of this function.
     */
    public int getInputDim() {
      return INPUT_DIM;
    }

    /**
     * Returns the output dimension of this function.
     */
    public int getOutputDim() {
      return OUTPUT_DIM;
    }

    /**
     * Evaluates this function at the supplied input.
     *
     * @param input a vector with {@link #getInputDim()} elements
     * @return a vector with {@link #getOutputDim()} elements
     */
    public double[] evaluate(double[] input) {
      double[] y = new double[OUTPUT_DIM];
      y[XPOS] = input[XPOS] + Math.cos(input[HEADING]) * input[TVEL];
      y[YPOS] = input[YPOS] + Math.sin(input[HEADING]) * input[TVEL];
      y[HEADING] = input[HEADING] + input[RVEL];
      y[TVEL] = input[TVEL_REL_NOISE] * tVelCtrl + input[TVEL_ABS_NOISE];
      y[RVEL] = input[RVEL_REL_NOISE] * rVelCtrl + input[RVEL_ABS_NOISE];
      return y;
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
      Matrix J = new Matrix(OUTPUT_DIM, INPUT_DIM);
      J.set(XPOS, XPOS, 1.0d);
      J.set(XPOS, HEADING, -Math.sin(input[HEADING]) * input[TVEL]);
      J.set(XPOS, TVEL, Math.cos(input[HEADING]));
      J.set(YPOS, YPOS, 1.0d);
      J.set(YPOS, HEADING, Math.cos(input[HEADING]) * input[TVEL]);
      J.set(YPOS, TVEL, Math.sin(input[HEADING]));
      J.set(HEADING, HEADING, 1.0d);
      J.set(HEADING, RVEL, 1.0d);
      J.set(TVEL, TVEL_REL_NOISE, tVelCtrl);
      J.set(TVEL, TVEL_ABS_NOISE, 1.0d);
      J.set(RVEL, RVEL_REL_NOISE, rVelCtrl);
      J.set(RVEL, RVEL_ABS_NOISE, 1.0d);
      return J.getArray();
    }
    
    /**
     * Returns the distribution over the noise variables for this motion
     * model.  
     */
    public Gaussian getNoiseModel() {
      return motionNoiseModel;
    }
  };

  /**
   * Returns the motion model.
   */
  public NoisyVectorFunction getMotionModel(double[] c) {
    return new MotionModel(c);
  }

  /**
   * A distribution over the noise variables of the odometry model.
   */
  protected Gaussian odometryNoiseModel;


  /**
   * Creates an odometry noise model.
   *
   * @param tVelRelNoise the variance of the relative noise in the
   *                     translational velocity odometry measurement
   * @param tVelAbsNoise the variance of the absolute noise in the
   *                     translational velocity odometry measurement
   *                     specified in (meters/second)<sup>2</sup>
   * @param rVelRelNoise the variance of the relative noise in the
   *                     rotational velocity odometry measurement
   * @param rVelAbsNoise the variance of the absolute noise in the
   *                     rotational velocity odometry measurement
   *                     specified in (radians/second)<sup>2</sup>
   * @return a Gaussian distribution over a vector variable whose
   *         components correspond to the relative and absolute errors 
   *         in the translational and rotational odometry measurements
   */
  protected static Gaussian getOdometryNoiseModel(double tVelRelNoise,
						  double tVelAbsNoise,
						  double rVelRelNoise,
						  double rVelAbsNoise) {
    Variable tVelRelNoiseVar = 
      new Variable("Relative noise in translational velocity odometry", 1);
    Variable tVelAbsNoiseVar = 
      new Variable("Absolute noise in translational velocity odometry", 1);
    Variable rVelRelNoiseVar = 
      new Variable("Relative noise in rotational velocity odometry", 1);
    Variable rVelAbsNoiseVar = 
      new Variable("Absolute noise in rotational velocity odometry", 1);
    ListSet noiseVars = new ListSet(tVelRelNoiseVar);
    noiseVars.add(tVelAbsNoiseVar);
    noiseVars.add(rVelRelNoiseVar);
    noiseVars.add(rVelAbsNoiseVar);
    Matrix mu = new Matrix(4, 1);
    mu.set(0, 0, 1.0d);
    mu.set(1, 0, 0.0d);
    mu.set(2, 0, 1.0d);
    mu.set(3, 0, 0.0d);
    Matrix sigma = new Matrix(4, 4);
    sigma.set(0, 0, tVelRelNoise);
    sigma.set(1, 1, tVelAbsNoise);
    sigma.set(2, 2, rVelRelNoise);
    sigma.set(3, 3, rVelAbsNoise);
    return new Gaussian(noiseVars, mu, sigma, true);
  }

  /**
   * Creates an odometry noise model using parameters specified by system
   * properties or default values (when the properties are not defined).
   *
   * @return a Gaussian distribution over a vector variable whose
   *         components correspond to the relative and absolute errors 
   *         in the translational and rotational odometry measurements
   */
  protected static Gaussian getOdometryNoiseModel() {
    double tVelRelNoise = 
	Double.parseDouble(System.getProperty("javaslam.ODO_TVELRELNOISE",
					      Double.toString(0.02 * 0.02)));
    double tVelAbsNoise = 
	Double.parseDouble(System.getProperty("javaslam.ODO_TVELABSNOISE",
					      Double.toString(0.01 * 0.01)));
    double rVelRelNoise = 
	Double.parseDouble(System.getProperty("javaslam.ODO_RVELRELNOISE",
					      Double.toString(0.04 * 0.04)));
    double rVelAbsNoise = 
	Double.parseDouble(System.
			   getProperty("javaslam.ODO_RVELABSNOISE",
				       Double.
				       toString(Math.pow(Math.PI / 360.0d, 2.0))));
    return getOdometryNoiseModel(tVelRelNoise,
				 tVelAbsNoise,
				 rVelRelNoise,
				 rVelAbsNoise);
  }

  /**
   * The odometry model of the robot.  The robot receives noisy
   * measurements of its translational and rotational velocity.  */
  public class OdometryModel 
    implements NoisyVectorFunction, ExtendedVectorFunction {

    /**
     * The input dimension of the odometry model.
     */
    protected final static int INPUT_DIM = 9;

    /**
     * The output dimension of the odometry model.
     */
    protected final static int OUTPUT_DIM = 2;
  
    /* These are indexes into the input vector that pull out specific
     * elements.
     */
    protected final static int XPOS = 0;
    protected final static int YPOS = 1;
    protected final static int HEADING = 2;
    protected final static int TVEL = 3;
    protected final static int RVEL = 4;
    protected final static int TVEL_REL_NOISE = 5;
    protected final static int TVEL_ABS_NOISE = 6;
    protected final static int RVEL_REL_NOISE = 7;
    protected final static int RVEL_ABS_NOISE = 8;

    /* These are indexes into the output vector that pull out specific
     * elements.
     */
    protected final static int ODO_TVEL = 0;
    protected final static int ODO_RVEL = 1;

    /**
     * Returns the input dimension of this function.
     */
    public int getInputDim() {
      return INPUT_DIM;
    }

    /**
     * Returns the output dimension of this function.
     */
    public int getOutputDim() {
      return OUTPUT_DIM;
    }

    /**
     * Evaluates this function at the supplied input.
     *
     * @param input a vector with {@link #getInputDim()} elements
     * @return a vector with {@link #getOutputDim()} elements
     */
    public double[] evaluate(double[] input) {
      double[] y = new double[OUTPUT_DIM];
      y[ODO_TVEL] = 
	input[TVEL_REL_NOISE] * input[TVEL] + input[TVEL_ABS_NOISE];
      y[ODO_RVEL] = 
	input[RVEL_REL_NOISE] * input[RVEL] + input[RVEL_ABS_NOISE];
      return y;
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
      Matrix J = new Matrix(OUTPUT_DIM, INPUT_DIM);
      J.set(ODO_TVEL, TVEL_REL_NOISE, input[TVEL]);
      J.set(ODO_TVEL, TVEL, input[TVEL_REL_NOISE]);
      J.set(ODO_TVEL, TVEL_ABS_NOISE, 1.0d);
      J.set(ODO_RVEL, RVEL_REL_NOISE, input[RVEL]);
      J.set(ODO_RVEL, RVEL, input[RVEL_REL_NOISE]);
      J.set(ODO_RVEL, RVEL_ABS_NOISE, 1.0d);
      return J.getArray();
    }
    
    /**
     * Returns the distribution over the noise variables for this motion
     * model.  
     */
    public Gaussian getNoiseModel() {
      return odometryNoiseModel;
    }
  };

  /**
   * The odometry model.  Since the model has no parameters, it is
   * instantiated once.
   */
  protected OdometryModel odometryModel;

  /**
   * Returns the odometry model.
   */
  public NoisyVectorFunction getOdometryModel() {
    return odometryModel;
  }

  /**
   * A distribution over the noise variables of the landmark
   * measurement model.  
   */
  protected Gaussian landmarkNoiseModel;

  /**
   * Creates a landmark measurement noise model.
   *
   * @param rangeAbsNoise the variance of the absolute noise in the
   *                      range measurement specified in 
   *                      (meters/second)<sup>2</sup>
   * @param rangeRelNoise the variance of the relative noise in the
   *                      range measurements
   * @param bearingNoise  the variance of the absolute noise in the
   *                      bearing measurements specified in
   *                      (radians/second)<sup>2</sup>
   * @return a Gaussian distribution over a vector variable whose
   *         components correspond to the relative and absolute errors 
   *         in the range and bearing measurements
   */
  protected static Gaussian getLandmarkNoiseModel(double rangeRelNoise,
						  double rangeAbsNoise,
						  double bearingNoise) {
    Variable rangeRelNoiseVar = 
      new Variable("Relative noise in range measurements", 1);
    Variable rangeAbsNoiseVar = 
      new Variable("Absolute noise in range measurements", 1);
    Variable bearingNoiseVar = 
      new Variable("Absolute noise in bearing measurements", 1);
    ListSet noiseVars = new ListSet(rangeRelNoiseVar);
    noiseVars.add(rangeAbsNoiseVar);
    noiseVars.add(bearingNoiseVar);
    Matrix mu = new Matrix(3, 1);
    mu.set(0, 0, 1.0d);
    mu.set(1, 0, 0.0d);
    mu.set(2, 0, 0.0d);
    Matrix sigma = new Matrix(3, 3);
    sigma.set(0, 0, rangeRelNoise);
    sigma.set(1, 1, rangeAbsNoise);
    sigma.set(2, 2, bearingNoise);
    return new Gaussian(noiseVars, mu, sigma, true);
  }

  /**
   * Creates a landmark measurement noise model using parameters
   * specified by system properties or default values (when the
   * properties are not defined).
   *
   * @return a Gaussian distribution over a vector variable whose
   *         components correspond to the relative and absolute errors 
   *         in the range and bearing measurements
   */
  protected static Gaussian getLandmarkNoiseModel() {
    double rangeRelNoise = 
	Double.parseDouble(System.getProperty("javaslam.RANGERELNOISE",
					      Double.toString(0.1 * 0.1)));
    double rangeAbsNoise = 
	Double.parseDouble(System.getProperty("javaslam.RANGEABSNOISE",
					      Double.toString(0.2 * 0.2)));
    double bearingNoise = 
	Double.parseDouble(System.getProperty("javaslam.BEARINGNOISE",
					      Double.
					      toString(Math.pow(Math.PI / 90.0d, 2.0))));
    return getLandmarkNoiseModel(rangeRelNoise,
				 rangeAbsNoise,
				 bearingNoise);
  }

  /**
   * The landmark measurement model of the robot.  The robot receives
   * noisy measurements of the bearing and range to landmarks.  */
  public class LandmarkMeasurementModel 
    implements NoisyVectorFunction, ExtendedVectorFunction {

    /**
     * The input dimension of the landmark measurement model.
     */
    protected final static int INPUT_DIM = 10;

    /**
     * The output dimension of the landmark measurement model.
     */
    protected final static int OUTPUT_DIM = 2;
  
    /* These are indexes into the input vector that pull out specific
     * elements.
     */
    protected final static int XPOS = 0;
    protected final static int YPOS = 1;
    protected final static int HEADING = 2;
    protected final static int TVEL = 3;
    protected final static int RVEL = 4;
    protected final static int LM_XPOS = 5;
    protected final static int LM_YPOS = 6;
    protected final static int RANGE_REL_NOISE = 7;
    protected final static int RANGE_ABS_NOISE = 8;
    protected final static int BEARING_NOISE = 9;

    /* These are indexes into the output vector that pull out specific
     * elements.
     */
    protected final static int DEL_X = 0;
    protected final static int DEL_Y = 1;

    /**
     * Returns the input dimension of this function.
     */
    public int getInputDim() {
      return INPUT_DIM;
    }

    /**
     * Returns the output dimension of this function.
     */
    public int getOutputDim() {
      return OUTPUT_DIM;
    }

    /**
     * Evaluates this function at the supplied input.
     *
     * @param input a vector with {@link #getInputDim()} elements
     * @return a vector with {@link #getOutputDim()} elements
     */
    public double[] evaluate(double[] input) {
      double range = 
	Math.sqrt(Math.pow(input[LM_XPOS] - input[XPOS], 2.0d) +
		  Math.pow(input[LM_YPOS] - input[YPOS], 2.0d));
      double noisyRange = input[RANGE_REL_NOISE] * range +
	input[RANGE_ABS_NOISE];
      double bearing = 
	Math.atan2(input[LM_YPOS] - input[YPOS],
		   input[LM_XPOS] - input[XPOS]) - input[HEADING];
      double noisyBearing = bearing + input[BEARING_NOISE];
      // Convert the noisy range-bearing measurements to cartesian
      // coordinates.
      double[] y = new double[OUTPUT_DIM];
      y[DEL_X] = noisyRange * Math.cos(noisyBearing);
      y[DEL_Y] = noisyRange * Math.sin(noisyBearing);
      return y;
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
      Matrix J = new Matrix(OUTPUT_DIM, INPUT_DIM);
      double range = 
	Math.sqrt(Math.pow(input[LM_XPOS] - input[XPOS], 2.0d) +
		  Math.pow(input[LM_YPOS] - input[YPOS], 2.0d));
      if (range == 0.0d) range = 1.0d;  // Avoid Inf
      double noisyRange = input[RANGE_REL_NOISE] * range +
	input[RANGE_ABS_NOISE];
      double bearing = 
	Math.atan2(input[LM_YPOS] - input[YPOS],
		   input[LM_XPOS] - input[XPOS]) - input[HEADING];
      double noisyBearing = bearing + input[BEARING_NOISE];
      // Calculate partial derivatives of the range.
      double dRANGE_dXPOS = -(input[LM_XPOS] - input[XPOS]) / range;
      double dRANGE_dLMXPOS = (input[LM_XPOS] - input[XPOS]) / range;
      double dRANGE_dYPOS = -(input[LM_YPOS] - input[YPOS]) / range;
      double dRANGE_dLMYPOS = (input[LM_YPOS] - input[YPOS]) / range;
      // Calculate partial derivatives of the noisy range.
      double dNRANGE_dXPOS = input[RANGE_REL_NOISE] * dRANGE_dXPOS;
      double dNRANGE_dLMXPOS = input[RANGE_REL_NOISE] * dRANGE_dLMXPOS;
      double dNRANGE_dYPOS = input[RANGE_REL_NOISE] * dRANGE_dYPOS;
      double dNRANGE_dLMYPOS = input[RANGE_REL_NOISE] * dRANGE_dLMYPOS;
      double dNRANGE_dRANGERELNOISE = range;
      double dNRANGE_dRANGEABSNOISE = 1.0d;
      // Calculate partial derivatives of the (noisy) bearing.
      double tmp = 1.0d / (1.0d + 
			   Math.pow((input[LM_YPOS] - input[YPOS]) / 
				    (input[LM_XPOS] - input[XPOS]), 2.0d));
      double dNBEARING_dLMYPOS = tmp / (input[LM_XPOS] - input[XPOS]);
      double dNBEARING_dYPOS = -dNBEARING_dLMYPOS;
      double dNBEARING_dLMXPOS = -dNBEARING_dLMYPOS * 
	(input[LM_YPOS] - input[YPOS]) / 
	Math.pow(input[LM_XPOS] - input[XPOS], 2.0d);
      double dNBEARING_dXPOS = -dNBEARING_dLMXPOS; 
      double dNBEARING_dHEADING = -1.0d;
      double dNBEARING_dBEARINGNOISE = 1.0d;
      // Calculate partial derivatives of the noisy X displacement.
      J.set(DEL_X, XPOS, dNRANGE_dXPOS * Math.cos(noisyBearing) 
	    - noisyRange * Math.sin(noisyBearing) * dNBEARING_dXPOS);
      J.set(DEL_X, LM_XPOS, dNRANGE_dLMXPOS * Math.cos(noisyBearing) 
	    - noisyRange * Math.sin(noisyBearing) * dNBEARING_dLMXPOS);
      J.set(DEL_X, YPOS, dNRANGE_dYPOS * Math.cos(noisyBearing) 
	    - noisyRange * Math.sin(noisyBearing) * dNBEARING_dYPOS);
      J.set(DEL_X, LM_YPOS, dNRANGE_dLMYPOS * Math.cos(noisyBearing) 
	    - noisyRange * Math.sin(noisyBearing) * dNBEARING_dLMYPOS);
      J.set(DEL_X, HEADING, 
	    -noisyRange * Math.sin(noisyBearing) * dNBEARING_dHEADING);
      J.set(DEL_X, RANGE_REL_NOISE, 
	    dNRANGE_dRANGERELNOISE * Math.cos(noisyBearing));
      J.set(DEL_X, RANGE_ABS_NOISE, 
	    dNRANGE_dRANGEABSNOISE * Math.cos(noisyBearing));
      J.set(DEL_X, BEARING_NOISE, 
	    -noisyRange * Math.sin(noisyBearing) * dNBEARING_dBEARINGNOISE);
      // Calculate partial derivatives of the noisy Y displacement.
      J.set(DEL_Y, XPOS, dNRANGE_dXPOS * Math.sin(noisyBearing) 
	    + noisyRange * Math.cos(noisyBearing) * dNBEARING_dXPOS);
      J.set(DEL_Y, LM_XPOS, dNRANGE_dLMXPOS * Math.sin(noisyBearing) 
	    + noisyRange * Math.cos(noisyBearing) * dNBEARING_dLMXPOS);
      J.set(DEL_Y, YPOS, dNRANGE_dYPOS * Math.sin(noisyBearing) 
	    + noisyRange * Math.cos(noisyBearing) * dNBEARING_dYPOS);
      J.set(DEL_Y, LM_YPOS, dNRANGE_dLMYPOS * Math.sin(noisyBearing) 
	    + noisyRange * Math.cos(noisyBearing) * dNBEARING_dLMYPOS);
      J.set(DEL_Y, HEADING, 
	    noisyRange * Math.cos(noisyBearing) * dNBEARING_dHEADING);
      J.set(DEL_Y, RANGE_REL_NOISE, 
	    dNRANGE_dRANGERELNOISE * Math.sin(noisyBearing));
      J.set(DEL_Y, RANGE_ABS_NOISE, 
	    dNRANGE_dRANGEABSNOISE * Math.sin(noisyBearing));
      J.set(DEL_Y, BEARING_NOISE, 
	    noisyRange * Math.cos(noisyBearing) * dNBEARING_dBEARINGNOISE);
      return J.getArray();
    }

    /**
     * Returns the distribution over the noise variables for this motion
     * model.  
     */
    public Gaussian getNoiseModel() {
      return landmarkNoiseModel;
    }
  };

  /**
   * The landmark measurement model.
   */
  public LandmarkMeasurementModel landmarkMeasurementModel;

  /**
   * Returns the landmark observation model.
   */
  public NoisyVectorFunction getMeasurementModel() {
    return landmarkMeasurementModel;
  }

  /**
   * The inverse landmark measurement model of the robot. 
   */
  public class InverseMeasurementModel 
    implements NoisyVectorFunction, ExtendedVectorFunction {

    /**
     * The input dimension of the model.
     */
    protected final static int INPUT_DIM = 8;

    /**
     * The output dimension of the model.
     */
    protected final static int OUTPUT_DIM = 2;
  
    /* These are indexes into the input vector that pull out specific
     * elements.
     */
    protected final static int XPOS = 0;
    protected final static int YPOS = 1;
    protected final static int HEADING = 2;
    protected final static int TVEL = 3;
    protected final static int RVEL = 4;
    protected final static int RANGE_REL_NOISE = 5;
    protected final static int RANGE_ABS_NOISE = 6;
    protected final static int BEARING_NOISE = 7;

    /* These are indexes into the output vector that pull out specific
     * elements.
     */
    protected final static int LM_XPOS = 0;
    protected final static int LM_YPOS = 1;

    /* These are indexes into the measurement vector that pull out specific
     * elements.
     */
    protected final static int DEL_X = 0;
    protected final static int DEL_Y = 1;

    /**
     * The noisy range measurement.
     */
    protected double range;

    /**
     * The noisy bearing measurement.
     */
    protected double bearing;

    /**
     * Constructor.
     *
     * @param z the measurement, represented as a displacement vector
     *          in the robot's cartesian coordinate frame; 
     *          <code>z[0]</code> is the frontal distance to the landmark, 
     *          and <code>z[0]</code> is the lateral distance to the landmark
     */
    public InverseMeasurementModel(double[] z) {
      range = Math.sqrt(z[DEL_X] * z[DEL_X] + z[DEL_Y] * z[DEL_Y]);
      bearing = Math.atan2(z[DEL_Y], z[DEL_X]);
    }
    
    /**
     * Returns the input dimension of this function.
     */
    public int getInputDim() {
      return INPUT_DIM;
    }

    /**
     * Returns the output dimension of this function.
     */
    public int getOutputDim() {
      return OUTPUT_DIM;
    }

    /**
     * Evaluates this function at the supplied input.
     *
     * @param input a vector with {@link #getInputDim()} elements
     * @return a vector with {@link #getOutputDim()} elements
     */
    public double[] evaluate(double[] input) {
      double[] y = new double[OUTPUT_DIM];
      y[LM_XPOS] = input[XPOS] + 
	(range / input[RANGE_REL_NOISE] - input[RANGE_ABS_NOISE]) *
	Math.cos(input[HEADING] + bearing - input[BEARING_NOISE]);
      y[LM_YPOS] = input[YPOS] + 
	(range / input[RANGE_REL_NOISE] - input[RANGE_ABS_NOISE]) *
	Math.sin(input[HEADING] + bearing - input[BEARING_NOISE]);
      return y;
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
      Matrix J = new Matrix(OUTPUT_DIM, INPUT_DIM);
      J.set(LM_XPOS, XPOS, 1.0d);
      J.set(LM_XPOS, RANGE_REL_NOISE,
	    (-range / Math.pow(input[RANGE_REL_NOISE], 2.0d)) *
	    Math.cos(input[HEADING] + bearing - input[BEARING_NOISE]));
      J.set(LM_XPOS, RANGE_ABS_NOISE,
	    -Math.cos(input[HEADING] + bearing - input[BEARING_NOISE]));
      J.set(LM_XPOS, HEADING, 
	    -(range / input[RANGE_REL_NOISE] - input[RANGE_ABS_NOISE]) *
	    Math.sin(input[HEADING] + bearing - input[BEARING_NOISE]));
      J.set(LM_XPOS, BEARING_NOISE, 
	    (range / input[RANGE_REL_NOISE] - input[RANGE_ABS_NOISE]) *
	    Math.sin(input[HEADING] + bearing - input[BEARING_NOISE]));
      J.set(LM_YPOS, YPOS, 1.0d);
      J.set(LM_YPOS, RANGE_REL_NOISE,
	    (-range / Math.pow(input[RANGE_REL_NOISE], 2.0d)) *
	    Math.sin(input[HEADING] + bearing - input[BEARING_NOISE]));
      J.set(LM_YPOS, RANGE_ABS_NOISE,
	    -Math.sin(input[HEADING] + bearing - input[BEARING_NOISE]));
      J.set(LM_YPOS, HEADING, 
	    (range / input[RANGE_REL_NOISE] - input[RANGE_ABS_NOISE]) *
	    Math.cos(input[HEADING] + bearing - input[BEARING_NOISE]));
      J.set(LM_YPOS, BEARING_NOISE, 
	    -(range / input[RANGE_REL_NOISE] - input[RANGE_ABS_NOISE]) *
	    Math.cos(input[HEADING] + bearing - input[BEARING_NOISE]));
      return J.getArray();
    }
    
    /**
     * Returns the distribution over the noise variables for this 
     * model.  
     */
    public Gaussian getNoiseModel() {
      return landmarkNoiseModel;
    }
  };

  /**
   * Returns the landmark observation model.
   */
  public NoisyVectorFunction getInverseMeasurementModel(double[] z) {
    return new InverseMeasurementModel(z);
  }
}
