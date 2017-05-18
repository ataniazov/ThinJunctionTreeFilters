package javaslam.tjt;

import java.util.*;
import javaslam.util.*;
import javaslam.prob.*;
import javaslam.tjt.graph.*;
import Jama.Matrix;

/**
 * A thin (approximate) junction tree for a Gaussian graphical model.
 * The width of this junction tree (and hence its time complexity for
 * inference) can be constrained by employing <i>variable
 * contractions</i>, in which variables are removed from clusters in
 * such a way as to preserve the running-intersection property.
 *
 * <p>In this implementation, the <i>width</i> of a junction tree is
 * defined to be the maximum <i>size</i> of its clusters; another
 * implementation could redefine the width to be the maximum
 * <i>dimension</i> of its clusters.  This latter notion of width
 * makes more sense because the complexity of inference is related
 * more tightly to the dimension of its clusters than their sizes;
 * however, it also requires splitting vector variables, which is
 * complex to implement.</p>
 * 
 * <p>This class records counts of all floating point operations using
 * {@link Flops#count(long)} (except those used in the service of
 * debugging and avoiding numerical errors).</p>
 *
 * @see "<i>Thin junction tree filters for simultaneous localization
 * and mapping.</i> Mark A. Paskin.  Computer Science Division
 * Technical Report CSD-02-1198, University of California, Berkeley.
 * 2002."
 * @author Mark A. Paskin
 * @version $Revision: 1.17 $ ($Date: 2003/01/10 07:22:05 $)
 */
public class ThinJunctionTree extends JunctionTree {
    
  /**
   * A priority queue of clusters prioritized by their size.
   */
  protected PriorityQueue clustersBySize;

  /**
   * Default constructor.
   */
  public ThinJunctionTree() {
    super();
    this.clustersBySize = new PriorityQueue();
  }

  /**
   * Returns a String representation of this thin junction tree.
   */
  public String toString() {
    return "thin " + super.toString();
  }

  /*------------------------------------------------------------------*
   * One of the functionalities added by this class over JunctionTree *
   * is to maintain a list of clusters sorted by their sizes (as *
   * a priority queue).  The methods below constitute the methods     *
   * for accessing and updating this list.                            *
   *------------------------------------------------------------------*/

  /**
   * Computes the width of this junction tree.  The width is the
   * size of the largest cluster in the junction tree.
   *
   * @return the size of the largest cluster in this junction tree
   */
  public int width() {
    if (clusters.isEmpty()) return 0;
    int width = ((Integer)clustersBySize.highest()).intValue();
    if (TJTF.blather) {
      TJTF.say("Width: " + width);
      TJTF.say("Priority queue: " + clustersBySize);
    }
    return width;
  }

  /**
   * Updates the {@link ThinJunctionTree#clustersBySize clustersBySize}
   * index of cluster sizes.
   *
   * @param c   the cluster whose size has changed
   */
  protected void updateClusterSize(Cluster c) {
    Integer s = new Integer(c.getSize());
    if (TJTF.blather) TJTF.say("Updating size of " + c + " to " + s);
    clustersBySize.reprioritize(c, s);
    if (TJTF.blather) TJTF.say("Priority queue: " + clustersBySize);
  }

  /**
   * Returns the cluster with largest size in the junction tree.
   *
   * @return the cluster with largest size in the junction tree
   */
  protected Cluster largestCluster() {
    Cluster largest = (Cluster)clustersBySize.peek();
    if (TJTF.debug && !clusters.contains(largest))
      throw new InternalError("Cluster queue contains orphan clusters ");
    return largest;
  }

  /*------------------------------------------------------------------*
   * So that the priority queue of clusters is always current, this   *
   * class extends some key methods of the JunctionTree class to      *
   * record updates in the priority queue.  These methods are below.  *
   *------------------------------------------------------------------*/

  /**
   * Extends the supplied cluster so that it contains the supplied
   * variable <i>without preserving consistency or the running-intersection 
   * property</i>.  This method is overrided here to call 
   * {@link JunctionTree#enlarge(JunctionTree.Cluster,Variable) enlarge}
   * and then record the change in the cluster's size.
   *
   * @param cluster a cluster of the junction tree to be enlarged
   * @param var     the variable to be added to the cluster
   */
  protected void enlarge(Cluster c, Variable var) {
    super.enlarge(c, var);
    updateClusterSize(c);
  }

  /**
   * Reduces the supplied cluster so that it no longer contains the
   * supplied variable.  The cluster potential is updated by
   * marginalizing out this variables as well.  This method preserves
   * consistency, but it makes no attempt to ensure the running
   * intersection property; that property will only persist if
   * <code>c</code> is a leaf of the subtree induced by
   * <code>var</code>.  This method is overrided here to call 
   * {@link JunctionTree#reduce(JunctionTree.Cluster,Variable) reduce}
   * and then record the change in the cluster's size.
   *
   * @param cluster a cluster of the junction tree to be enlarged
   * @param var     the variable to be added to the cluster
   */
  protected void reduce(Cluster c, Variable var) {
    super.reduce(c, var);
    updateClusterSize(c);
  }

  /**
   * Creates an empty cluster and attaches it as a leaf off of the
   * supplied cluster.  This method is overrided here to call 
   * {@link JunctionTree#newLeaf(JunctionTree.Cluster) newLeaf}
   * and then record the cluster in the priority queue.
   *
   * @param c the cluster to which the new cluster will be attached as
   *          a leaf; this may be <code>null</code> if the junction tree
   *          currently has no clusters
   * @return the new cluster
   * @throws IllegalArgumentException if <code>c == null</code>
   *                                  but this junction tree is not empty
   */
  protected Cluster newLeaf(Cluster c) {
    Cluster newCluster = super.newLeaf(c);
    updateClusterSize(newCluster);
    return newCluster;
  }

  /**
   * Removes a cluster from this junction tree.  This method only
   * works if <code>c</code> is not connected to any other clusters.
   * This method is overrided here to call {@link
   * JunctionTree#remove(JunctionTree.Cluster) remove} and then remove
   * the cluster from the priority queue.
   *
   * @param c the cluster to be removed
   * @throws IllegalArgumentException if <code>c</code> is connected
   *                                  to any other clusters of the 
   *                                  junction tree
   */
  protected void remove(Cluster c) {
    super.remove(c);
    clustersBySize.remove(c);
  }

  /*-----------------------------------------------------------------*
   * The following methods are general methods required to keep a    *
   * junction tree thin: cluster cloning, variable contraction,      *
   * contraction error estimation, etc.                              *
   *-----------------------------------------------------------------*/

  /**
   * Represents a <i>variable contraction</i> of a variable <i>v</i>
   * along an edge <code>e</code> = <i>C<sub>i</sub> ->
   * C<sub>j</sub></i>.  If this contraction is executed, then the
   * <i>v</i> is marginalized out of the cluster <i>C<sub>j</sub></i>
   * and also from the separator potential of <code>e</code>.  If
   * there are other edges emanating from <i>C<sub>j</sub></i> that
   * carry <i>v</i>, then <i>v</i> is first (recursively) contracted
   * along these edges.  This operation preserves the singly-connected
   * and running-intersection properties (as well as consistency), but
   * violates the potential property.  Thus, it results in a junction
   * tree that is valid for some graphical model, but not that of the
   * original junction tree.
   */
  public class Contraction {
    /**
     * The variable to contracted.
     */
    public Variable v;

    /**
     * The edge along which the variable is to be contracted.
     */
    public JTEdge e;

    /**
     * If not <code>null</code>, this is a list of {@link
     * ThinJunctionTree.Contraction}s that must first be executed
     * before this contraction can be executed.
     */
    protected List descendants;

    /**
     * The approximation error induced by this contraction.
     */
    protected double err;

    /**
     * Initializer.
     *
     * @param e     an edge in the junction tree 
     * @param v     a variable contained in the cluster to which 
     *              <code>e</code> is directed
     * @param local if <code>true</code>, then this contraction assumes that
     *              the cluster to which <code>e</code> is directed is a 
     *              leaf in the subtree induced by <code>v</code> and does
     *              not perform any recursive contractions
     * @throws IllegalArgumentException if either of the clusters incident 
     *                                  to <code>e</code> do not contain 
     *                                  <code>v</code>
     */
    protected void initialize(JTEdge e, Variable v, boolean local) {
      if (TJTF.blather) TJTF.say("Initializing " + (local ? "local " : "") + 
		       "edge contraction of " + v + " along " + e);
      Cluster c = (Cluster)e.to;
      Cluster d = (Cluster)e.from;
      /*
      if (d.getVariables().containsAll(c.getVariables()))
	throw new InternalError("Cluster subsumed by contracting neighbor");
      */
      if (!c.contains(v) || !d.contains(v))
	throw new IllegalArgumentException("Cannot contract variable along " +
					   "edge when indicent clusters do " + 
					   "not contain it");
      this.e = e;
      this.v = v;
      this.descendants = new LinkedList();
      if (!local) {
	// Add all descendant edge contractions.
	for (Iterator i = c.edges().iterator(); i.hasNext();) {
	  JTEdge f = (JTEdge)i.next();
	  if (f.getReverse().equals(e)) continue;
	  if (((Cluster)f.to).contains(v)) {
	    if (TJTF.blather) TJTF.say("Adding descendant edge contraction of " + 
			     v + " along " + f);
	    descendants.add(new Contraction(f, v));
	  }
	}
      }
    }

    /**
     * Constructor.
     *
     * @param e     an edge in the junction tree 
     * @param v     a variable contained in the cluster to which 
     *              <code>e</code> is directed
     * @param local if <code>true</code>, then this contraction assumes that
     *              the cluster to which <code>e</code> is directed is a 
     *              leaf in the subtree induced by <code>v</code> and does
     *              not perform any recursive contractions
     * @throws IllegalArgumentException if either of the clusters incident 
     *                                  to <code>e</code> do not contain 
     *                                  <code>v</code>
     */
    public Contraction(JTEdge e, Variable v, boolean local) {
      initialize(e, v, local);
      this.err = computeError();
    }

    /**
     * Constructor.
     *
     * @param e an edge in the junction tree 
     * @param v a variable contained in the cluster to which 
     *          <code>e</code> is directed
     * @throws IllegalArgumentException if either of the clusters incident 
     *                                  to <code>e</code> do not contain 
     *                                  <code>v</code>
     */
    public Contraction(JTEdge e, Variable v) {
      this(e, v, false);
    }

    /**
     * Constructor.
     *
     * @param e     an edge in the junction tree 
     * @param v     a variable contained in the cluster to which 
     *              <code>e</code> is directed
     * @param local if <code>true</code>, then this contraction assumes that
     *              the cluster to which <code>e</code> is directed is a 
     *              leaf in the subtree induced by <code>v</code> and does
     *              not perform any recursive contractions
     * @param error the error of this contraction
     * @throws IllegalArgumentException if either of the clusters incident 
     *                                  to <code>e</code> do not contain 
     *                                  <code>v</code>
     */
    public Contraction(JTEdge e, Variable v, 
		       boolean local, double error) {
      initialize(e, v, local);
      this.err = error;
    }

    /**
     * Returns the approximation error that will result if this
     * contraction is executed, which is the relative entropy (or
     * Kullback-Lieblier divergence) from the original distribution to
     * the resulting distribution.
     *
     * @return the approximation error that will result if this
     * contraction is executed
     */
    public double error() {
      return this.err;
    }

    /**
     * Executes the edge contraction.
     */
    public void execute() {
      if (TJTF.debug) checkValid();
      if (TJTF.verbose) TJTF.say("Executing edge contraction of " + v + " along " + e);
      // Recursively contract along the descendant edges.
      for (Iterator i = descendants.iterator(); i.hasNext();) {
	Contraction ec = (Contraction)i.next();
	ec.execute();
      }
      /* Marginalize the variable out of the terminal cluster and the
       * separator potential. */
      ThinJunctionTree.this.reduce((Cluster)e.to, v);
      e.getPotential().marginalizeOut(Collections.singleton(v));
      Cluster c = (Cluster)e.to;
      if (false) {
	/* If the cluster from which the variable was contracted is
	 * subsumed by one of its neighbors, then merge its subsumer
	 * into it. (This means it will still be around for a recursive
	 * contraction.) */
	Cluster d = c.getSubsumingNeighbor();
	if (d != null) {
	  if (TJTF.debug && d.equals((Cluster)e.from))
	    throw new InternalError("Cluster subsumed by contracting neighbor");
	  merge((JTEdge)c.getEdge(d));
	}
      }
      if (TJTF.debug) checkValid();
      if (TJTF.debug && !clusters.contains(c))
	throw new InternalError("Cluster orphaned by edge contraction");
    }

    /**
     * Computes the approximation error introduced by a variable
     * contraction along an edge <code>e</code> = <i>C<sub>i</sub></i>
     * -> <i>C<sub>j</sub></i>.  This approximation error is the
     * Kullback-Liebler divergence from the original junction tree's
     * distribution to the distribution that would result if {@link
     * #execute()} were called.
     *
     * <p>If the terminus <i>C<sub>j</sub></i> is a leaf of the
     * subtree induced by <code>v</code>, then the approximation error
     * is the conditional mutual information
     * <nobr><b>I</b>(<i>C<sub>j</sub></i> - <i>S</i>; <code>v</code>
     * | <i>S</i>)</nobr> where <i>S</i> is the separator of
     * <code>e</code>.  If <i>C<sub>j</sub></i> is not a leaf of the
     * subtree induced by <code>v</code>, then the approximation error
     * is the sum of the contraction errors of <code>v</code> along
     * all edges <i>C<sub>j</sub></i> -> <i>C<sub>k</sub></i> (except
     * <i>C<sub>j</sub></i> -> <i>C<sub>i</sub></i>).</p>
     *
     * @return      the Kullback-Liebler divergence from the original 
     *              junction tree's distribution to the distribution
     *              that would result if {@link #execute()} were called
     */
    public double computeError() {
      ThinJunctionTree.this.passFlows();
      Cluster c = (Cluster)e.to;
      double sumError = 0.0d;
      for (Iterator i = descendants.iterator(); i.hasNext();) {
	Contraction ec = (Contraction)i.next();
	sumError += ec.error();
      }
      Gaussian cp = c.getPotential();
      Gaussian sp = e.getPotential();
      ListSet ls = new ListSet(v);
      Matrix lcx = cp.getLambda(ls, null);
      Matrix lsx = sp.getLambda(ls, null);
      Flops.count(Flops.det(lcx.getRowDimension(), true) +
		  Flops.det(lsx.getRowDimension(), true) +
		  2 * Flops.log() + 3);
      return sumError + (0.5d) * (Math.log(lcx.det()) - Math.log(lsx.det()));
    }

    /**
     * Returns <code>true</code> if this and the supplied object are
     * equal.  Two edge contractions are equal if their edges and
     * variables are the same.
     *
     * @return <code>true</code> if this and the supplied object are
     *         equal
     */
    public boolean equals(Object o) {
      return ((o instanceof Contraction) &&
	      ((Contraction)o).e.equals(e) &&
	      ((Contraction)o).v.equals(v));
    }

    /**
     * Returns the hash code of this edge contraction.
     *
     * @return the hash code of this edge contraction
     */
    public int hashCode() {
      return e.hashCode() + v.hashCode();
    }
  }

  /**
   * Creates a variable contraction of the supplied variable along the
   * supplied edge.
   *
   * @param v the variable to be contracted
   * @param e the edge along which the variable should be contracted
   * @return a variable contraction of <code>v</code> along <code>e</code>
   */
  public Contraction getContraction(JTEdge e, Variable v) {
    return new Contraction(e, v, false);
  }

  /**
   * Duplicates the supplied cluster and attaches the duplicate as a
   * leaf off of the original cluster.  This operation preserves
   * consistency and validity.
   *
   * @param c the cluster to be cloned
   * @return the clone
   */
  public Cluster clone(Cluster c) {
    if (TJTF.verbose) TJTF.say("Cloning " + c);
    Cluster clone = newLeaf(c);
    JTEdge e = (JTEdge)clone.getEdge(c);
    /* We have to use 'enlarge' here rather than 'extend' since the
     * latter choice would create the clone, and then remove it
     * because it is nonmaximal!
     */
    for (Iterator i = c.getVariables().iterator(); i.hasNext();)
	ThinJunctionTree.this.enlarge(clone, (Variable)i.next());
    clone.getPotential().set(c.getPotential());
    e.getPotential().set(c.getPotential());
    return clone;
  }


  /*-----------------------------------------------------------------*
   * The following classes and methods implement the thinning        *
   * strategy.                                                       *
   *-----------------------------------------------------------------*/

  /**
   * Thins the junction tree via a sequence of variable contractions
   * so that no clusters have a size greater than <code>limit</code>.
   *
   * @param limit the cluster size limit
   * @return the approximation error introduced by the thinning; this
   *         error is the Kullback-Liebler divergence from the
   *         original distribution to the thinner distribution (in
   *         natural logarithmic units)
   */
  protected double thin(int limit) {
    double error = 0.0;
    if (TJTF.blather) TJTF.say("Thinning...");
    while (width() > limit) {
      Cluster largest = largestCluster();
      error += thin(largest, limit, null);
    }
    if (TJTF.blather) TJTF.say("done thinning.");
    return error;
  }
    
  /**
   * Determines if the supplied cluster is a non-root leaf of the
   * subtree induced by the supplied variable.  The subtree induced by
   * <code>v</code> is the set of clusters that contain
   * <code>v</code>; because of the singly-connected and
   * running-intersection properties, this set must be a tree.  If
   * <code>c</code> is a leaf of this subtree (and is not the only
   * cluster in this subtree), then this method returns the {@link
   * JunctionTree.JTEdge JTEdge} to <code>c</code> from its neighbor
   * in the subtree; otherwise, this method returns <code>null</code>.
   *
   * <p><b>Algorithm:</b> This implementation tests if <code>c</code>
   * has exactly one neighbor that also contains <code>v</code>.
   * Assuming the singly-connected and running-intersection properties
   * hold, this is necessary and sufficient for <code>c</code> to be a
   * leaf of <code>v</code>'s induced subtree.  Then it returns the
   * edge to <code>c</code> from this neighbor.</p>
   *
   * @param c a cluster
   * @param v a variable contained in <code>c</code>
   * @return the {@link JunctionTree.JTEdge JTEdge} to
   *         <code>c</code> from its neighbor in the subtree, or
   *         <code>null</code> if <code>c</code> is not a non-root leaf
   *         of the subtree induced by <code>v</code>
   * @throws IllegalArgumentException if <code>c</code> does not
   *                                  contain <code>v</code>
   */
  protected JTEdge isSubtreeLeaf(Cluster c, Variable v) {
    if (!c.contains(v)) 
      throw new IllegalArgumentException("" + c + " does not contain " + v);
    JTEdge edge = null;
    for (Iterator i = c.edges().iterator(); i.hasNext();) {
      JTEdge e = (JTEdge)i.next();
      if (((Cluster)e.to).contains(v)) {
	if (edge != null)
	  return null;
	else
	  edge = e.getReverse();
      }
    }
    if (TJTF.blather) 
      TJTF.say(c.toString() + " is a leaf of the subtree induced by " + v);
    return edge;
  }

  /**
   * Removes variables from the supplied cluster so as to ensure its
   * size is no greater than <code>limit</code>.  The approximation
   * error introduced is returned.  If not enough variables can be
   * removed to satisfy the supplied limit, then an exception is thrown.
   *
   * <p>This method evaluates the approximation error that would arise
   * from contracting each variable that does not only occur in
   * <code>c</code> and chooses the contraction that results in the
   * least amount of error.</p>
   *
   * @param c       the cluster to be thinned
   * @param limit   the size to which the cluster should be thinned
   * @param exclude if not <code>null</code>, this variable will not
   *                be removed from <code>c</code>
   * @return the approximation error introduced by the thinning; this
   *         error is the Kullback-Liebler divergence from the
   *         original distribution to the thinner distribution (in
   *         natural logarithmic units) 
   * @throws UnsupportedOperationException if not enough variables in
   *         <code>c</code> appear elsewhere
   */
  public double thin(Cluster c, int limit, Variable exclude) {
    if (TJTF.debug) checkValid();
    if (TJTF.debug && !clusters.contains(c))
      throw new InternalError("Attempt to thin an orphaned cluster");
    if (TJTF.verbose) TJTF.say("Thinning " + c);
    double err = 0.0d;
    while (c.getSize() > limit) {
      /* Each time we perform a contraction from this cluster, it is
       * necessary to recompute the approximation errors that would be
       * induced by all other possible contractions.  This is because
       * the membership of the cluster has changed (and possibly those
       * of distant clusters as well) , and therefore the mutual
       * informations will be different. */
      Contraction best = null;
      for (Iterator i = c.getVariables().iterator(); i.hasNext();) {
	Variable v = (Variable)i.next();
	if (v.equals(exclude)) continue;
	if (getClustersWith(v).size() > 1) {
	  /* This variable occurs elsewhere.  If this cluster is a
	   * leaf of the subtree induced by v, then examine whether
	   * contracting v from this cluster is better than any
	   * contractions found thus far.
	   */
	  JTEdge e = isSubtreeLeaf(c, v);
	  if (e != null) {
	    Contraction other = new Contraction(e, v);
	    if ((best == null) || (other.err < best.error()))
	      best = other;
	    continue;
	  }
	  /* This variable occurs elsewhere and this cluster is an
	   * internal node of the subtree induced by v; thus, there
	   * are at least two edges along which we could contract v.
	   * We examine each to see if it is better than the best
	   * contraction found thus far.  We avoid redundant
	   * computation in computing the errors of each possibility
	   * as follows.  We first compute the error of contracting v
	   * to c independently from each neighbor carrying v.  Then,
	   * to compute the cost of contracting v from c along e = (c,
	   * d), we sum the errors of contracting v to c along all
	   * v-carrying neighbors (except d) with the local cost of
	   * contracting along e.
	   */
	  HashSet contractions = new HashSet();
	  double sumError = 0.0;
	  for (Iterator j = c.edges().iterator(); j.hasNext();) {
	    e = (JTEdge)j.next();
	    Cluster d = (Cluster)e.to;
	    if (!d.contains(v)) continue;
	    Contraction ec = new Contraction(e, v);
	    sumError += ec.error();
	    contractions.add(ec);
	  }
	  if (TJTF.debug && !clusters.contains(c))
	    throw new InternalError("Attempt to thin an orphaned cluster");
	  // There must be at least two contractions.
	  if (contractions.size() < 2)
	    throw new InternalError("Cluster index or edges inconsistent");
	  // Now compute the error of contracting v from c along each
	  // of the edges from neighbors carrying v.
	  for (Iterator j = contractions.iterator(); j.hasNext();) {
	    Contraction ec = (Contraction)j.next();
	    e = ec.e.getReverse();
	    double error = sumError - ec.error() + 
	      new Contraction(e, v, true).error();
	    if ((best == null) || (error < best.error()))
	      best = new Contraction(e, v, false, error);
	  }
	}
      }
      // Every cluster supports one contraction as long as its size is
      // at least three.
      if (best == null)
	throw new UnsupportedOperationException("Cluster supports " + 
						"no contractions.");
      // Now perform the contraction.
      if (TJTF.debug && !clusters.contains(c))
	throw new InternalError("Attempt to thin an orphaned cluster");
      if (TJTF.debug) checkValid();
      best.execute();
      if (TJTF.debug) checkValid();
      if (TJTF.debug && !clusters.contains(c))
	  throw new InternalError("Cluster orphaned by contraction: " + best);
      err += best.error();
    }
    return err;
  }

  /**
   * Contracts <code>var</code> from all clusters but <code>cluster</code>
   * (which must contain <code>var</code>).
   *
   * @param var the variable to be contracted
   * @param cluster a cluster containing var
   */
  public void contractTo(Variable var, Cluster cluster) {
    if (!cluster.contains(var))
      throw new IllegalArgumentException("Cluster does not contain variable");
    for (Iterator j = cluster.edges().iterator(); j.hasNext();) {
      JTEdge e = (JTEdge)j.next();
      Cluster d = (Cluster)e.to;
      if (!d.contains(var)) continue;
      Contraction ec = new Contraction(e, var, false, Double.NaN);
      ec.execute();
    }
  }

  /**
   * Contracts <code>var</code> from all clusters but one and returns
   * that cluster.
   *
   * @param var the variable to be contracted
   * @return the only cluster containing var 
   */
  public Cluster contract(Variable var) {
    /* First we construct a queue of (edges incident to) leaves of
     * the subtree induced by var. */
    List unexamined = new LinkedList();
    JTEdge e = null;
    for (Iterator i = getClustersWith(var).iterator(); i.hasNext();) 
      if ((e = isSubtreeLeaf((Cluster)i.next(), var)) != null) 
	unexamined.add(e);
    // Now we repeatedly process the leaves until there are no more.
    PriorityQueue contractions = new PriorityQueue();
    while (getClustersWith(var).size() > 1) {
      // Process all of the unexamined leaves. 
      while (!unexamined.isEmpty()) {
	e = (JTEdge)unexamined.remove(0);
	Contraction ec = new Contraction(e, var);
	contractions.enqueue(ec, new Double(-ec.error()));
      }
      /* Perform the edge contraction that induces the smallest amount
       * of error, and then add that leaf's neighbor to the queue of
       * leaves to be processed. */
      Contraction ec = (Contraction)contractions.dequeue();
      ec.execute();
      Cluster neighbor = (Cluster)ec.e.from;
      if ((e = isSubtreeLeaf(neighbor, var)) != null) 
	unexamined.add(e);
    }
    return (Cluster)getClustersWith(var).iterator().next();
  }
}
