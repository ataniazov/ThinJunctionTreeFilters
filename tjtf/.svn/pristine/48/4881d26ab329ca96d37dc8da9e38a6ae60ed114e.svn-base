package javaslam.slam;

import java.util.*;
import javaslam.prob.*;

/**
 * <p>An interface implemented by filters for the Simultaneous
 * Localization and Mapping (<font size="-1">SLAM</font>) problem.
 * Such classes must assign landmarks integer identifiers and support
 * queries on the states of the robot and landmarks.</p>
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.10 $ ($Date: 2003/02/07 02:26:18 $) */
public interface SLAMFilter {

  /**
   * Gets the state variable associated with the robot.
   */
  public Variable getRobotVariable();

  /**
   * Gets the state variable associated with the landmark with the
   * supplied identifier.
   */
  public Variable getLandmarkVariable(int id);
  
  /**
   * Gets the identifier associated with the supplied landmark state
   * variable.
   */
  public int getLandmarkId(Variable lm);

  /**
   * Returns the number of landmarks known to this filter.
   */
  public int getNumLandmarks();

  /**
   * Adds a new landmark to the map and returns its state variable.
   */
  public Variable addLandmark(int id, int dim);

  /**
   * Returns the IDs of landmarks known to this filter.
   */
  public int[] getLandmarkIds();

  /**
   * Return <code>true</code> if this filter contains the landmark
   * with the supplied identifier.
   */
  public boolean contains(int id);

  /**
   * Returns the filtered estimate of the robot's state.
   *
   * @return the filtered estimate of the robot's state
   */
  public double[] getRobotEstimate();

  /**
   * Returns the filtered estimate of a landmark's state.
   *
   * @param id the identifier of the landmark
   * @return the filtered estimate of the landmark's state
   */
  public double[] getLandmarkEstimate(int id);
}
