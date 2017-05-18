package javaslam.slam;

import java.util.*;
import Jama.*;
import javaslam.util.*;
import javaslam.prob.*;
import javaslam.filter.*;
import javaslam.tjt.*;
import javaslam.tjt.graph.*;

/**
 * A thin junction tree filter for the Simultaneous Localization and
 * Mapping (SLAM) problem.
 *
 * @see "<i>Thin Junction Tree Filters for Simultaneous Localization
 *      and Mapping</i>. Mark Paskin.  Technical Report
 *      UCB/CSD-02-1198, University of California, Berkeley, 2002."
 * @author Mark A. Paskin
 * @version $Revision: 1.10 $ ($Date: 2003/02/07 20:06:31 $) 
 */
public class TJTSLAMFilter extends LGSLAMFilter {

  /**
   * The thin junction tree used to represent the filtered belief
   * state.
   */
  protected ThinJunctionTree jt;

  /**
   * The cluster size limit of the thin junction tree.  If this is
   * <code>-1</code>, then the junction tree is not thinned.
   */
  protected int width;

  /**
   * When the active cluster is cloned to admit a new landmark
   * variable, this is the size of its separator.
   */
  protected int overlap;

  /**
   * For convenience, a list-set containing only {@link #x}.
   */
  private final ListSet xSet;

  /**
   * Constructor.
   *
   * @param width   the maximum number of variables permitted in a
   *                single cluster of the thin junction tree; this 
   *                must be two or greater
   * @param overlap the size of separators for clones; this must be
   *                positive and less than width
   * @param mu      a <i>n</i>-vector containing the starting state 
   *                of the robot
   * @param sigma   an <i>n</i>-by-<i>n</i> positive semidefinite 
   *                covariance matrix giving the uncertainty in the 
   *                robot's initial state
   */
  public TJTSLAMFilter(int width,
		       int overlap,
		       double[] mu, 
		       double[][] sigma) {
    super(mu.length);
    jt = new ThinJunctionTree();
    xSet = new ListSet();
    xSet.add(x);
    Gaussian p = new Gaussian(xSet, mu, sigma, true);
    jt.times(p.reparameterize(true));
    if (width < 2)
      throw new IllegalArgumentException("Width must be two or greater");
    if ((overlap <= 0) || (overlap >= width))
      throw new IllegalArgumentException("Overlap must be positive " + 
					 "and less than width");
    this.width = width;
    this.overlap = overlap;
  }

  /**
   * Gets the junction tree representation of the belief state.
   */
  public JunctionTree getJunctionTree() {
    return jt;
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
   * covariance <code>R</code>.  Given the actual measurement
   * <code>z</code>, this method updates the belief state as follows:
   *
   * <ol type="1">
   *
   *    <li>A measurement potential is constructed over the state
   *    variable of the robot state and that of the observed
   *    landmark.</li>
   *
   *    <li>If the landmark has not yet been observed and the active
   *    cluster is full (meaning its size is the width limit), then:
   *    <ol type="a">
   *      <li>the active cluster is cloned;</li>
   *      <li>the clone is made the active cluster</li>
   *      <li>the robot state variable is contracted to the active
   *      cluster</li>
   *      <li>the active cluster is thinned until its separator
   *      contains <code>overlap</code> variables</li>
   *    </ol>
   *    Otherwise, the cluster closest to the active cluster that
   *    contains the state variable of the observed landmark is made
   *    the active cluster.</li>
   *
   *    <li>The measurement potential is multiplied into the active
   *    cluster.</li>
   *
   * </ol>
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
    Gaussian p = InformationFilter.
      getConditionalGaussian(vars, 
			     new Matrix(z0, z0.length), 
			     CD,
			     new Matrix(R), 
			     new Matrix(z, z.length));
    // Multiply in the potential.
    JunctionTree.Cluster best = null;
    if (jt.contains(lm)) {
      if (TJTF.verbose) 
	TJTF.say("Landmark " + lm + " previously observed.");
      // Select the best cluster as the cluster containing lm that
      // is closest to a cluster containing x.
      Set clusters = jt.getClustersWith(lm);
      int shortest = Integer.MAX_VALUE;
      for (Iterator i = clusters.iterator(); i.hasNext();) {
	JunctionTree.Cluster c = (JunctionTree.Cluster)i.next();
	if (c.contains(x)) {
	  best = c;
	  break;
	}
	// Find the path to the nearest cluster containing 'var'.
	List path = c.getShortestPath(new NodeFilter(x) {
	    public boolean satisfies(Node n) {
	      return ((JunctionTree.Cluster)n).contains((Variable)data);}});
	if (path == null)
	  throw new InternalError("Cannot find path to " + x);
	if ((best == null) || (path.size() < shortest)) {
	  best = c;
	  shortest = path.size();
	}
      }
      if (TJTF.verbose) 
	TJTF.say("Selected best cluster as " + best);
    } else {
      if (TJTF.verbose) 
	TJTF.say("Landmark " + lm + " newly observed.");
      // Find the smallest cluster containing x.
      Set clusters = jt.getClustersWith(x);
      int smallest = Integer.MAX_VALUE;
      for (Iterator i = clusters.iterator(); i.hasNext();) {
	JunctionTree.Cluster c = (JunctionTree.Cluster)i.next();
	if ((best == null) || (c.getSize() < smallest)) {
	  best = c;
	  smallest = c.getSize();
	  break;
	}
      }
      // If the best cluster cannot admit a new landmark variable
      // without exceeding the size limit, then clone it and make the
      // clone the best cluster.
      if ((width != -1) && (best.getSize() == width)) {
	jt.passFlows();
	JunctionTree.Cluster clone = jt.clone(best);
	jt.contractTo(x, clone);
	jt.thin(clone, overlap, x);
	best = clone;
      }
    }
    // Multiply the measurement potential into the best cluster.
    jt.times(p, best);
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
   * <code>Q</code>.  This method updates the belief state as follows:
   *
   * <ol type="1">
   *
   *    <li>An evolution potential is constructed over the robot state
   *    variables at the current and the next time step.</li>
   *
   *    <li>This robot's state variable at the current time step is
   *    contracted to the active cluster.</li>
   *
   *    <li>The evolution potential is multiplied into the active
   *    cluster.</li>
   *
   *    <li>This robot's state variable at the current time step is
   *    marginalized out of the active cluster.</li>
   *
   *    <li>The robot's state variable at the next time step is
   *    renamed to be the robot's state variable at the current time
   *    step.</li>
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
    // Create a variable that represents x at the next time step.
    Variable xt = new Variable(x.label + "'", x.dim);
    ListSet nextVars = new ListSet();
    nextVars.add(xt);
    // Construct the potential.
    Gaussian p = InformationFilter.
      getUnconditionalGaussian(nextVars, xSet, 
			       new Matrix(x0, x0.length), 
			       new Matrix(A), new Matrix(Q));
    // Contract the robot's state variable to a single cluster.
    jt.passFlows();
    JunctionTree.Cluster c = jt.contract(x);
    // Multiply in the evolution potential.
    jt.times(p, c);
    // Perform 'roll-up'.
    jt.marginalizeOut(x);
    jt.rename(xt, x);
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
   * <code>S</code>.  An odometry measurement potential over <i>x</i>
   * is constructed and multiplied into the active cluster.
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
    Gaussian p = InformationFilter.
      getConditionalGaussian(xSet, 
			     new Matrix(y0, y0.length), 
			     new Matrix(B), 
			     new Matrix(S), 
			     new Matrix(y, y.length));
    jt.times(p);
  }

  /**
   * Returns the filtered marginal potential (in the moment
   * parameterization) over the robot's state.
   *
   * @return the filtered marginal potential over the robot's state.
   */
  public Gaussian getRobotMarginal() {
    Gaussian p = jt.getMarginal(xSet, true).reparameterize(true);
    return p;
  }

  /**
   * Returns the filtered marginal potential over a landmark's state
   * (in the moment parameterization).
   *
   * @param ids the landmark identifier
   * @return the filtered marginal potential over a landmark's state
   */
  public Gaussian getLandmarkMarginal(int id) {
    return jt.getMarginal(Collections.singleton(getLandmarkVariable(id)),
			  true).reparameterize(true);
  }

  /**
   * This method computes the joint filtered distribution over the
   * states of the robot and a landmark.  The distribution may be a
   * product-of-marginals approximation if it is cheaper to compute.
   *
   * @param id the landmark's identifier
   * @return the joint filtered distribution over the
   *         states of the robot and a landmark (in the moment 
   *         parameterization)
   */
  public Gaussian getRobotLandmarkMarginal(int id) {
    ListSet tmp = new ListSet(getRobotVariable());
    tmp.add(getLandmarkVariable(id));
    return jt.getMarginal(tmp, false).reparameterize(true);
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
    Set vars = new HashSet();
    for (int i = 0; i < ids.length; i++)
      vars.add(getLandmarkVariable(ids[i]));
    Map var2marg = jt.getMarginals(vars);
    Gaussian[] p = new Gaussian[ids.length];
    for (int i = 0; i < ids.length; i++)
      p[i] = (Gaussian)var2marg.get(getLandmarkVariable(ids[i]));
    return p;
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
    Set vars = new HashSet();
    for (int i = 0; i < ids.length; i++)
      vars.add(getLandmarkVariable(ids[i]));
    HashMap map = new HashMap();
    // Get as many true marginals as possible.
    Set clusters = jt.getClustersWith(x);
    for (Iterator i = clusters.iterator(); i.hasNext();) {
      JunctionTree.Cluster c = (JunctionTree.Cluster)i.next();
      HashSet intersection = new HashSet(c.getVariables());
      intersection.retainAll(vars);
      if (!intersection.isEmpty()) {
	int u = c.getPotential().getDimension();
	Gaussian p = c.getPotential().reparameterize(false);
	ListSet tmp = new ListSet(getRobotVariable());
	for (Iterator j = intersection.iterator(); j.hasNext();) {
	  Variable landmark = (Variable)j.next();
	  tmp.add(landmark);
	  map.put(landmark, p.marginalize(tmp, false));
	  tmp.remove(landmark);
	  vars.remove(landmark);
	}
      }
    }
    // Compute the rest using a product of marginals approximation.
    if (!vars.isEmpty()) {
      Gaussian p = 
	jt.getMarginal(Collections.singleton(x), true).reparameterize(true);
      Map rest = jt.getMarginals(vars);
      for (Iterator i = vars.iterator(); i.hasNext();) {
	Variable lm = (Variable)i.next();
	Gaussian q = (Gaussian)rest.get(lm);
	map.put(lm, q.times(p, true));
      }
    }
    Gaussian[] p = new Gaussian[ids.length];
    for (int i = 0; i < ids.length; i++)
      p[i] = (Gaussian)map.get(getLandmarkVariable(ids[i]));
    return p;
  }
}
