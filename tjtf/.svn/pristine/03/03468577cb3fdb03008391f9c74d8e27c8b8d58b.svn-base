package javaslam.slam;

import java.util.*;
import javaslam.prob.*;

/**
 * <p>An abstract class extended by classes that implement filters for
 * the Simultaneous Localization and Mapping (SLAM) problem.  This
 * class implements a simple landmark indexing scheme in which each
 * landmark is represented by an integer identifier.</p>
 *
 * <p>A more sophisticated implementation of this class could
 * implement a cache of landmark marginals.</p>
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.1 $ ($Date: 2003/02/07 02:26:18 $) 
 */
public abstract class AbstractSLAMFilter implements SLAMFilter {

  /**
   * The variable <i>x</i>, which represents the state of the robot.
   */
  public final Variable x;

  /**
   * A map from {@link Integer} landmark identifiers to {@link
   * Variable}s representing the corresponding landmark's state.
   */
  protected Map id2lm;

  /**
   * A map from {@link Variable}s representing landmarks' states
   * to {@link Integer} identifiers of the corresponding landmarks .
   */
  protected Map lm2id;

  /**
   * Constructor.
   *
   * @param dim the dimension of the robot state variable
   */
  public AbstractSLAMFilter(int dim) {
    x = new Variable("X", dim);
    id2lm = new HashMap();
    lm2id = new HashMap();
  }

  /**
   * Gets the state variable associated with the robot.
   */
  public final Variable getRobotVariable() {
    return x;
  }

  /**
   * Gets the state variable associated with the landmark with the
   * supplied identifier.
   */
  public Variable getLandmarkVariable(int id) {
    return (Variable)id2lm.get(new Integer(id));
  }
  
  /**
   * Gets the identifier associated with the supplied landmark state
   * variable.
   */
  public int getLandmarkId(Variable lm) {
    return ((Integer)lm2id.get(lm)).intValue();
  }

  /**
   * Returns the number of landmarks known to this filter.
   */
  public int getNumLandmarks() {
    return id2lm.size();
  }

  /**
   * Adds a new landmark to the map and returns its state variable.
   */
  public Variable addLandmark(int id, int dim) {
    Variable lm = new Variable("Landmark #" + id, dim);
    Integer idInt = new Integer(id);
    id2lm.put(idInt, lm);
    lm2id.put(lm, idInt);
    return lm;
  }

  /**
   * Returns the IDs of landmarks known to this filter.
   */
  public int[] getLandmarkIds() {
    int[] ids = new int[id2lm.size()];
    int k = 0;
    for (Iterator i = id2lm.keySet().iterator(); i.hasNext();)
      ids[k++] = ((Integer)i.next()).intValue();
    return ids;
  }

  /**
   * Return <code>true</code> if this filter contains the landmark
   * with the supplied identifier.
   */
  public boolean contains(int id) {
    return id2lm.containsKey(new Integer(id));
  }


  /**
   * Returns the filtered estimate of the robot's state.
   *
   * @return the filtered estimate of the robot's state
   */
  public abstract double[] getRobotEstimate();

  /**
   * Returns the filtered estimate of a landmark's state.
   *
   * @param id the identifier of the landmark
   * @return the filtered estimate of the landmark's state
   */
  public abstract double[] getLandmarkEstimate(int id);
}
