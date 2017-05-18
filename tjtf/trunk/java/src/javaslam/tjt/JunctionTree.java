package javaslam.tjt;

import java.util.*;
import javaslam.prob.*;
import javaslam.util.*;
import javaslam.tjt.graph.*;

/**
 * A junction tree for a Gaussian graphical model.  A junction tree
 * over a set of potentials is a graph whose nodes (called
 * <i>clusters</i>) are sets of variables and that has three
 * properties:
 *
 * <ol>
 *   <li><i>singly-connected property</i>: the graph is a tree; and,</li>
 *   <li><i>potential property</i>: all potentials that have been 
 *       added to this junction tree are covered by at least one 
 *       cluster; and,</li>
 *   <li><i>running intersection property</i>: if <i>C<sub>1</sub></i> and
 *       <i>C<sub>2</sub></i> are clusters that contain the variable
 *       <i>x</i>, then all clusters (and separators) on the (unique) path
 *       between <i>C<sub>1</sub></i> and <i>C<sub>2</sub></i> also 
 *       contain <i>x</i>.</li>
 * </ol>
 *
 * In addition, each edge has an associated <i>separator</i> which is
 * the intersection of its two incident clusters.  Finally, every
 * cluster and separator has an associated potential function over its
 * variables.  A junction tree is <i>consistent</i> if each separator
 * potential agrees with (the appropriate marginal of) its incident
 * clusters' potentials.  If a junction tree is <i>valid</i> (meaning
 * all of the above properties hold) and it is consistent, then each
 * of its cluster (and separator) potentials are marginal
 * distributions for the graphical model defined by the normalized
 * product of the potentials in the junction tree.
 *
 * <p>This implementation is designed to efficiently support
 * incremental operations such as adding new variables, multiplying in
 * new potentials, and marginalizing out variables; at all times, the
 * junction tree remains valid and consistent so that cluster
 * marginals are accessible in constant time.</p>
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.17 $ ($Date: 2002/12/28 21:59:22 $)
 */
public class JunctionTree {
    
  /**
   * The {@link Variable Variable}s present in this junction tree.
   * This set's iterator is consistent with the order in which the
   * variables were added to the junction tree.
   */
  protected ListSet variables;

  /**
   * The {@link Cluster Cluster}s in this junction tree.
   */
  protected ListSet clusters;

  /**
   * Maps each {@link Variable} to the {@link Set Set} of {@link
   * Cluster Cluster}s containing it.  Variables that are not
   * currently in the junction tree may still be in this map, but they
   * will map to empty sets.
   */
  protected Map varToClusters;

  /**
   * The set of clusters that have been updated with new evidence but
   * have not yet distributed their evidence (due to lazy or adaptive
   * message passing).  When this set is empty, the junction tree is
   * consistent.
   */
  protected Set updated;

  /**
   * A threshold used in adaptive message passing.  If
   * <code>significance</code> is non-negative, then when distributing
   * evidence from a cluster, messages are only propagated while the
   * differential relative entropy (or {@link Gaussian#kl(Gaussian)
   * kl} divergence) from a separator potential to its new value (in
   * nats) after the message is passed is larger than
   * <code>significance</code>.
   *
   * @see JunctionTree.Cluster#distributeEvidence(JunctionTree.Cluster,double)
   */
  protected double significance;

  /**
   * Default constructor.  <code>significance</code> is initialized to
   * <code>-1.0</code> so that all messages will be passed.
   */
  public JunctionTree() {
    clusters = new ListSet();
    variables = new ListSet();
    varToClusters = new HashMap();
    updated = new HashSet();
    significance = -1.0d;
  }

  /**
   * Sets a threshold used in adaptive message passing.  If
   * <code>s</code> is non-negative, then when distributing evidence
   * from a cluster, messages are only propagated while the
   * differential relative entropy (or {@link Gaussian#kl(Gaussian)
   * kl} divergence) from a separator potential to its
   * new value (in nats) after the message is passed is larger than
   * <code>s</code>.
   *
   * @param s the new threshold value
   * @see JunctionTree.Cluster#distributeEvidence(JunctionTree.Cluster,double) 
   */
  public void setSignificance(double s) {
    significance = s;
  }

  /**
   * Returns <code>true</code> if this junction tree contains the
   * supplied variable.
   */
  public boolean contains(Variable var) {
    return variables.contains(var);
  }

  /**
   * Gets an unmodifiable set of the {@link Variable Variable}s in
   * this junction tree.  The iteration order of this set is the order
   * in which the variables were added to the juntion tree.
   */
  public Set getVariables() {
    return Collections.unmodifiableSet(variables);
  }

  /**
   * Returns a {@link java.util.Set}s of the {@link Cluster}s of this
   * junction tree.
   *
   * @return a {@link java.util.Set}s of the {@link Cluster}s
   *         of this junction tree.
   */
  public Set getClusters() {
    return Collections.unmodifiableSet(clusters);
  }

  /**
   * Removes a variable from this junction tree.
   *
   * @param var the variable to be removed
   * @throws IllegalArgumentException if <code>var</code> is present in any 
   *                                  cluster of the junction tree
   */
  protected void remove(Variable var) {
    if (!getClustersWith(var).isEmpty())
      throw new IllegalArgumentException("Cannot remove " + var + 
					 " while it is present in clusters");
    varToClusters.remove(var);
    variables.remove(var);
  }
      
  /**
   * Adds a variable to this junction tree.
   *
   * @param var the variable to be added
   */
  protected void add(Variable var) {
    if (variables.contains(var)) return;
    variables.add(var);
    varToClusters.put(var, new HashSet());
  }

  /**
   * Renames a variable in this junction tree.
   *
   * @param var   the original variable
   * @param subst the variable substituted for the original variable
   * @throws IllegalArgumentException if <code>var</code> and
   *         <code>subst</code> have differing dimensions
   */
  public void rename(Variable var, Variable subst) {
    if (!variables.contains(var)) return;
    Set s = (Set)varToClusters.get(var);
    for (Iterator i = s.iterator(); i.hasNext();)
      ((Cluster)i.next()).rename(var, subst);
    varToClusters.remove(var);
    varToClusters.put(subst, s);
    variables.remove(var);
    variables.add(subst);
  }

  /**
   * Returns an unmodifiable set of {@link Cluster Cluster}s
   * containing the supplied variable.
   *
   * @return the set of {@link Cluster Cluster}s of this junction tree
   *         that contain <code>var</code> or <code>null</code> if 
   *         <code>var</code> is not in the junction tree
   */
  public Set getClustersWith(Variable var) {
    Set s = (Set)varToClusters.get(var);
    if (s != null)
      return Collections.unmodifiableSet(s);
    else 
      return null;
  }

  /**
   * Returns the smallest cluster containing the supplied variables, or 
   * <code>null</code> if there is no such cluster.
   *
   * @param vars a set of {@link Variable Variable}s
   * @return the smallest cluster containing <code>vars</code> or 
   *         <code>null</code> if there is no such cluster
   */
  protected Cluster getCover(Set vars) {
    if (!variables.containsAll(vars)) return null;
    Cluster c = bestCover(vars);
    if (c.contains(vars)) return c;
    else return null;
  }

  /**
   * Multiplies a new potential into a particular cluster of this
   * junction tree and restores validity (but not consistency).  If
   * the supplied potential acts on variables that are not yet present
   * in the junction tree, these variables are first added to the
   * junction tree.
   *
   * @param p       the potential to be multiplied in 
   * @param cluster the cluster into which it should be multiplied
   */
  public void times(Gaussian p, Cluster cluster) {
    if (cluster.getJunctionTree() != this) 
      throw new IllegalArgumentException("Cluster does not " + 
					 "belong to this junction tree.");
    // Make sure the cluster contains all the potential's variables.
    for (Iterator i = p.getVariables().iterator(); i.hasNext();)
      extend(cluster, (Variable)i.next());
    if (TJTF.verbose) TJTF.say("Multiplying potential over " + 
		     p.getVariables() + " into " + cluster);
    cluster.getPotential().times(p, true);
    // Record that this cluster needs to distribute evidence.
    updated.add(cluster);
    // Check to see if this cluster now subsumes any of its neighbors.
    Cluster n = cluster.getSubsumedNeighbor();
    if (n != null) {
      if (TJTF.verbose) TJTF.say("Multiplying potential over " + 
		       p.getVariables() + " into " + cluster + 
		       " caused it to subsume " + n);
      JTEdge f = (JTEdge)cluster.getEdge(n);
      merge(f);
    }
  }

  /**
   * Multiplies a new potential into the junction tree and restores
   * validity (but not consistency).  If the supplied potential acts
   * on variables that are not yet present in the junction tree, these
   * variables are first added to the junction tree.
   *
   * @param p the potential to be multiplied in
   * @return the cluster into which <code>p</code> was multiplied
   */
  public Cluster times(Gaussian p) {
    // Find or create a cluster that covers the variables of p.
    Cluster cluster = getCover(p.getVariables());
    if (cluster == null) {
      if (TJTF.verbose) TJTF.say("Creating cover for potential over " + 
		       p.getVariables());
      cluster = createCover(p.getVariables());
    } else
      if (TJTF.verbose) TJTF.say("Found cover " + cluster + 
	  " for potential over " + p.getVariables());
    // Multiply in the potential.
    times(p, cluster);
    return cluster;
  }

  /**
   * Extends the supplied cluster so that it contains the supplied
   * variable <i>without preserving consistency or the
   * running-intersection property</i>.  
   *
   * <p><i>It is part of the contract of this class that whenever a
   * variable is added to a cluster, this method must be
   * invoked.</i></p>
   *
   * @param cluster a cluster of the junction tree to be enlarged
   * @param var     the variable to be added to the cluster
   * @see JunctionTree#extend(JunctionTree.Cluster,Variable)
   */
  protected void enlarge(Cluster c, Variable var) {
    if (c.contains(var)) return;
    if (TJTF.blather) TJTF.say("Enlarging " + c + " to include " + var);
    add(var);
    c.extend(Collections.singleton(var));
    ((Set)varToClusters.get(var)).add(c);
  }
    
  /**
   * Reduces the supplied cluster so that it no longer contains the
   * supplied variable.  The cluster potential is updated by
   * marginalizing out this variables as well.  This method preserves
   * consistency, but it makes no attempt to ensure the running
   * intersection property; that property will only persist if
   * <code>c</code> is a leaf of the subtree induced by
   * <code>var</code>.
   *
   * <p><i>It is part of the contract of this class that whenever a
   * variable is removed from a cluster, this method must be invoked;
   * However, this method is not invoked if a cluster is removed from
   * the junction tree entirely; in that case, {@link
   * JunctionTree#remove(JunctionTree.Cluster) remove} is called
   * instead.</i></p>
   *
   * @param cluster a cluster of the junction tree to be enlarged
   * @param var     the variable to be added to the cluster
   */
  protected void reduce(Cluster c, Variable var) {
    c.marginalizeOut(Collections.singleton(var));
    ((Set)varToClusters.get(var)).remove(c);
  }

  /**
   * Removes a cluster from this junction tree.  This method only
   * works if <code>c</code> is not connected to any other clusters.
   *
   * <p><i>It is part of the contract of this class that this method
   * is invoked whenever a cluster is removed from the junction
   * tree.</i></p>
   *
   * @param c the cluster to be removed
   * @throws IllegalArgumentException if <code>c</code> is connected
   *                                  to any other clusters of the 
   *                                  junction tree
   */
  protected void remove(Cluster c) {
    if (!c.edges().isEmpty())
      throw new IllegalArgumentException("Cannot remove cluster with " + 
					 "neighbors");
    for (Iterator i = c.getVariables().iterator(); i.hasNext();)
      ((Set)varToClusters.get((Variable)i.next())).remove(c);
    clusters.remove(c);
  }

  /**
   * Creates an empty cluster and attaches it as a leaf off of the
   * supplied cluster.  This method is the <i>only</i> correct way to
   * create a new junction tree cluster.</p>
   *
   * <p><i>It is part of the contract of this class that this method
   * is invoked whenever a new cluster is added to the junction
   * tree.</i></p>
   *
   * @param c the cluster to which the new cluster will be attached as
   *          a leaf; this may be <code>null</code> if the junction tree
   *          currently has no clusters
   * @return the new cluster
   * @throws IllegalArgumentException if <code>c == null</code>
   *                                  but this junction tree is not empty
   */
  protected Cluster newLeaf(Cluster c) {
    Cluster newCluster = new Cluster();
    if (c != null)
      link(c, newCluster);
    else if (!clusters.isEmpty())
      throw new IllegalArgumentException("Must attach new cluster");
    clusters.add(newCluster);
    return newCluster;
  }

  /**
   * <p>Minimally alters the structure and parameterization of the
   * junction tree so that <code>cluster</code> covers
   * <code>var</code> and validity and consistency are preserved.</p>
   *
   * <p><b>Algorithm:</b> In order to restore the running intersection
   * property, the closest cluster <i>C</i> containing
   * <code>var</code> is found, and <code>var</code> is added to all
   * clusters (and separators) on the path to <i>C</i>.  Consistency
   * is restored by passing flows <i>backwards</i> along this path.
   * This operation can result in non-maximal clusters; these are
   * subsequently detected by traversing the path and removed by
   * cluster merging.</p>
   *
   * @param cluster a cluster of the junction tree to be extended
   * @param var     the variable to be added to the cluster
   */
  protected void extend(Cluster cluster, Variable var) {
    if (cluster.contains(var)) return;
    // Check if the variable currently exists in the junction tree.
    if (variables.contains(var)) {
      // Find the path to the nearest cluster containing 'var'.
      List path = cluster.getShortestPath(new NodeFilter(var) {
	  public boolean satisfies(Node n) {
	    return ((Cluster)n).contains((Variable)data);}});
      // Since the variable exists in the model, there must be a path.
      if (path == null)
	throw new InternalError("Cannot find path to variable " + var);
      if (TJTF.verbose) TJTF.say("Shortest path to another cluster containing " + 
		       var + ": " + path);
      /* To restore the running intersection property, add 'var' to
       * all clusters on the path. 
       */
      Collections.reverse(path);
      for (Iterator j = path.iterator(); j.hasNext();) {
	JTEdge e = (JTEdge)j.next();
	// Extend the cluster.
	enlarge((Cluster)e.from, var);
	// Extend the potential of the separator.
	e.getPotential().extend(Collections.singleton(var));
	/* Pass a flow along the reverse of this edge; if 'from' was
	 * previously consistent, then this will restore its consistency. */
	if (TJTF.verbose) TJTF.say("Restoring consistency of " + e.from);
	e.getReverse().passFlow();
	if (TJTF.debug && !e.consistent())
	  throw new InternalError("Consistency not restored!");
      }
      /* The process of restoring the running intersection property
       * above can lead to non-maximal clusters.  To eliminate these
       * clusters, we examine each updated cluster C on the path to
       * see if it subsumes one of its neighbors.  If it does, then we
       * merge them together.
       */
      if (TJTF.verbose) TJTF.say("Searching for nonmaximal clusters along " + path);
      for (Iterator j = path.iterator(); j.hasNext();) {
        JTEdge e = (JTEdge)j.next();
	Cluster d = (Cluster)e.from;
	if (!clusters.contains(d)) continue;
	Cluster n = d.getSubsumedNeighbor();
	if (n != null) {
	  /* If this is not the last cluster on the reverse path
	   * (which is the cluster we are extending), then we merge
	   * the subsumer into the subsumee.  (This ensures that all
	   * clusters on the remainder of the path still exist.) If
	   * this is 'cluster', however, we merge the subsumed into
	   * the subsumer.  This guarantees to the caller that
	   * 'cluster' is still in the junction tree.
	   */
	  if (j.hasNext()) {
	    JTEdge f = (JTEdge)n.getEdge(d);
	    merge(f);
	  } else {
	    JTEdge f = (JTEdge)d.getEdge(n);
	    merge(f);
	  }
	}
      }
    }
    // Now that the running intersection property will hold, extend
    // the cluster.
    enlarge(cluster, var);
  }
  
  /**
   * Returns the cluster whose intersection with the supplied set of
   * variables is largest and whose size is the smallest.
   *
   * @param vars a {@link java.util.Set} of {@link Variable}s
   * @return the cluster whose intersection with <code>vars</code> is
   *         largest and whose size is the smallest
   */
  protected Cluster bestCover(Set vars) {
    /* Compute the set of clusters whose intersection with 'vars' is
     * nonempty. */
    HashSet possibilities = new HashSet();
    for (Iterator i = vars.iterator(); i.hasNext();) {
      Variable v = (Variable)i.next();
      Set s = getClustersWith(v);
      if (s != null) possibilities.addAll(s);
    }
    /* Of the possibilities, find that set whose intersection with
     * 'vars' is the largest (and of those, select the one whose
     * size is the smallest). */
    Cluster best = null;
    int bestOverlap = 0;
    for (Iterator i = possibilities.iterator(); i.hasNext();) {
      Cluster d = (Cluster)i.next();
      Set s = d.getVariables();
      int overlap = 0;
      for (Iterator j = s.iterator(); j.hasNext();)
	overlap += vars.contains(j.next()) ? 1 : 0;
      if ((best == null) || 
	  (overlap > bestOverlap) ||
	  ((overlap == bestOverlap) && (d.getSize() < best.getSize()))) {
	best = d;
	bestOverlap = overlap;
      }
    }
    return best;
  }

  /**
   * Returns the cluster that, if {@link
   * #extend(JunctionTree.Cluster,Variable) extend}ed with the
   * variables in <code>vars</code>, would cause the fewest number of
   * cluster {@link #enlarge(JunctionTree.Cluster,Variable)
   * enlarge}ments.  This differs from the result of {@link
   * #bestCover(Set)} because it takes into account the cost of
   * restoring the running intersection property.  (Note that there
   * may be a better cluster to extend---in that its extension will
   * increase the width of the junction tree by a smaller amount---but
   * this criterion is a fairly good and cheap heuristic.)
   *
   * @param vars a {@link java.util.Set} of {@link Variable}s
   * @return the best choice of clusters to {@link
   *         #extend(JunctionTree.Cluster,Variable) extend} with the 
   *         variables in <code>vars</code>
   */
  protected Cluster bestCoverToExtend(Set vars) {
    /* Compute the set of clusters whose intersection with 'vars' is
     * nonempty. */
    HashSet possibilities = new HashSet();
    for (Iterator i = vars.iterator(); i.hasNext();) {
      Variable var = (Variable)i.next();
      Set s = getClustersWith(var);
      if (s != null) possibilities.addAll(s);
    }
    if (possibilities.isEmpty())
      throw new InternalError("No clusters contain any of " + vars);
    /* Of the possibilities, find that cluster that results in the 
     * smallest number of cluster enlargements.
     */
    Cluster best = null;
    int bestNumEnlargements = 0;
    for (Iterator i = possibilities.iterator(); i.hasNext();) {
      Cluster c = (Cluster)i.next();
      int numEnlargements = 0;
      for (Iterator j = vars.iterator(); j.hasNext();) {
	Variable var = (Variable)j.next();
	if (!variables.contains(var)) continue;
	if (c.contains(var)) continue;
	// Find the path to the nearest cluster containing 'var'.
	List path = c.getShortestPath(new NodeFilter(var) {
	    public boolean satisfies(Node n) {
	      return ((Cluster)n).contains((Variable)data);}});
	numEnlargements += path.size();
      }
      if (TJTF.blather)
	TJTF.say("Extending " + c + " with " + vars + 
	    " requires " + numEnlargements + " enlargements.");
      if ((best == null) || (numEnlargements < bestNumEnlargements)) {
	best = c;
	bestNumEnlargements = numEnlargements;
      }
    }
    return best;
  }

  /**
   * Updates this junction tree so that it has a cluster containing
   * the supplied set of variables, while preserving validity and
   * consistency.  The {@link #bestCoverToExtend(Set)
   * bestCoverToExtend} is found and then is {@link
   * JunctionTree#extend(JunctionTree.Cluster,Variable) extend}ed to
   * cover each variable in <code>vars</code>.
   *
   * @param vars the set of {@link Variable Variable}s that must be covered 
   * @return the cover created for <code>vars</code>
   */
  protected Cluster createCover(Set vars) {
    // Choose the cluster to extend.
    Cluster c = null;
    if (clusters.isEmpty()) 
      c = newLeaf(null);
    else {
      Cluster best = bestCoverToExtend(vars);
      if (TJTF.verbose) TJTF.say("Best cluster to extend for " + vars + ": " + best);
      /* How to extend the best cluster:
       *   * If all variables of 'vars' are in the junction tree 
       *     already, then we enlarge 'best' to accomodate 'vars'.  (If
       *     we attached a new leaf off of 'best', then restoring the
       *     running intersection property would make 'best' and its new
       *     leaf identical.)
       *   * If 'vars' contains new variables, then:
       *     - if 'vars' contains 'best', then add 'vars' to 'best'
       *     - otherwise, 'vars' contains existing variables that 
       *       do not yet cohabitate; add a new leaf off of 'best' and 
       *       enlarge it to accomodate 'vars'.
       */
      if (variables.containsAll(vars)) {
	if (TJTF.verbose) TJTF.say("All the variables are already in the model.");
	// There must be at least one cluster, and so best != null
	c = best;
      } else if (vars.containsAll(best.getVariables())) {
	if (TJTF.verbose) TJTF.say("New variables are being added to the model " + 
	    "without new couplings of current variables");
	// There must be at least one cluster, and so best != null
	c = best;
      } else {
	if (TJTF.verbose) TJTF.say("New variables are being added to the model " + 
			 "with new couplings of current variables");
	// Create a new leaf off of best.
	c = newLeaf(best);
	if (TJTF.verbose) TJTF.say("Created a new cluster for " + vars + 
			 " as a leaf off of " + best);
      }
    }
    // First, enlarge the cluster to contain all the new variables.
    for (Iterator i = vars.iterator(); i.hasNext();) {
      Variable var = (Variable)i.next();
      if (!variables.contains(var)) enlarge(c, var);
    }
    // Finally, enlarge the cluster to contain all the old variables.
    for (Iterator i = vars.iterator(); i.hasNext();) {
      Variable var = (Variable)i.next();
      if (vars.contains(var)) extend(c, var);
    }
    if (TJTF.debug) checkValid();
    return c;
  } 

  /**
   * Merges the two clusters that are incident to the supplied edge.
   * The "to" cluster is merged into the "from" cluster.  All edges
   * incident to the "to" cluster are "swung" over to the "from"
   * cluster.  This method preserves validity and consistency.
   *
   * @param e the edge whose incident clusters are to be merged
   * @return the merged cluster (which is the same as <code>e.from</code>)
   */
  protected Cluster merge(JTEdge e) {
    Cluster from = (Cluster)e.from;
    Cluster to = (Cluster)e.to;
    if (TJTF.verbose) TJTF.say("Merging " + to + " into " + from);
    // Extend the 'from' cluster.
    for (Iterator i = to.getVariables().iterator(); i.hasNext();)
      enlarge(from, (Variable)i.next());
    // Update the potential from the 'from' cluster to become the
    // potential of the merged cluster (the product of the cluster
    // potentials, divided by the separator potential).
    if (TJTF.blather) TJTF.say("Updating the potential over " + from);
    Gaussian p = from.getPotential();
    p.times(to.getPotential(), true);
    p.div(e.getPotential(), true);
    // Swing all the edges from the "to" cluster over to the "from" cluster.
    if (TJTF.blather) TJTF.say("Swinging edges from " + to + " to " + from);
    to.swingEdgesTo(from);
    // Remove the 'to' cluster from the junction tree.
    if (TJTF.blather) TJTF.say("Removing " + to + " from the junction tree.");
    remove(to);
    return from;
  }

  /**
   * Merges all clusters containing a particular variable.  The
   * resulting merged cluster is returned.
   *
   * @param var the variable whose cover clusters are to be merged
   * @return the merged cluster
   */
  protected Cluster mergeClustersWith(Variable var) {
    Set s = getClustersWith(var);
    // This is implemented via a sequence of pairwise merges.
    if (s.isEmpty()) return null;
    while (s.size() > 1) {
      // Find an edge between two clusters containing var
      Cluster c = (Cluster)s.iterator().next();
      JTEdge edge = null;
      for (Iterator i = c.edges().iterator(); i.hasNext();) {
	JTEdge e = (JTEdge)i.next();
	if (((Cluster)e.to).contains(var)) {
	  edge = e;
	  break;
	}
      }
      if (edge == null) 
	throw new InternalError("Running intersection property violated");
      // Merge the two clusters incident to this edge.
      merge(edge);
    }
    return (Cluster)s.iterator().next();
  }

  /**
   * Marginalizes a variable out of this junction tree.  If the
   * junction tree was previously consistent, then it will be
   * consistent after this method is invoked.  All clusters containing
   * the variable are merged together and then the variable is
   * marginalized out of the resulting cluster.
   *
   * @param var the variable to be marginalized out
   */
  public void marginalizeOut(Variable var) {
    if (TJTF.debug && !consistent()) 
      throw new InternalError("Inconsistent before marginalizing out " + var);

    // Merge all clusters containing the variable.
    Cluster c = mergeClustersWith(var);
    // Now marginalize the variable out of the merged cluster.
    reduce(c, var);
    // Remove the variable from the set of this junction tree's variables.
    variables.remove(var);
    // Remove the variable from the cluster index.
    varToClusters.remove(var);

    // If the cluster is now subsumed by one of its neighbors, merge them.
    Cluster n = c.getSubsumingNeighbor();
    if (n != null) {
      JTEdge e = (JTEdge)c.getEdge(n);
      merge(e);
    }

    if (TJTF.debug && !consistent())
      throw new InternalError("Inconsistent after marginalizing out " + var);
  }

  /**
   * Restores consistency; all clusters into new potentials have been
   * multiplied that have not yet distributed their evidence do so.
   * If no new potentials have been multiplied into the junction tree
   * since the last time this method was invoked, then this method
   * does nothing.
   */
  public void passFlows() {
    while (!updated.isEmpty()) {
      Cluster c = (Cluster)updated.iterator().next();
      if (significance >= 0.0)
	c.distributeEvidence(null, significance);
      else
	c.distributeEvidence(null);
      updated.remove(c);
    }
  }

  /**
   * Extracts the marginal from the junction tree.  This method
   * first calls {@link #passFlows()} to ensure the junction tree is
   * consistent.
   *
   * @param vars  the set of {@link Variable Variable}s whose marginal 
   *              is to be computed
   * @param force if <code>true</code>, then the junction tree is
   *              restructured if necessary to compute the marginal; 
   *              otherwise, <code>null</code> is returned if there 
   *              is no cover for <code>vars</code> in this junction tree.
   * @return a marginal potential over <code>vars</code> (in the
   *         canonical parameterization), or <code>null</code> if
   *         <code>force</code> is <code>false</code> and
   *         <code>vars</code> do not reside together in a cluster of
   *         the junction tree
   */
  public Gaussian getMarginal(Set vars, boolean force) {
    passFlows();
    Cluster c = getCover(vars);
    if (c == null) {
      if (force)
	c = createCover(vars);
      else
	return null;
    }
    return c.getPotential().marginalize(vars, false);
  }

  /**
   * Extracts a set of unary marginals from the junction tree without
   * inverting any cluster potential more than once.  This method
   * first calls {@link #passFlows()} to ensure the junction tree is
   * consistent.
   *
   * @param  vars a collection of {@link Variable Variable}s, or 
   *         <code>null</code> to indicate all variables in this
   *         junction tree
   * @return a map whose keys are the (distinct) elements of
   *         <code>vars</code> and whose values are the corresponding 
   *         {@link Gaussian} marginals (in the moment parameterization)
   */
  public Map getMarginals(Collection vars) {
    passFlows();
    if (vars == null) vars = variables;
    HashMap marginals = new HashMap(vars.size());
    HashSet remaining = new HashSet(vars);
    if (!variables.containsAll(vars))
      throw new IllegalArgumentException("Not all variables in the model");
    while (!remaining.isEmpty()) {
      Cluster c = bestCover(remaining);
      if (TJTF.blather) TJTF.say("Best cover for " + remaining + ": " + c);
      Gaussian p = c.getPotential().reparameterize(false);
      for (Iterator i = p.getVariables().iterator(); i.hasNext();) {
	Variable v = (Variable)i.next();
	if (remaining.contains(v)) {
	  marginals.put(v, p.marginalize(Collections.singleton(v), false));
	  remaining.remove(v);
	}
      }
    }
    return marginals;
  }

  /**
   * Returns a String representation of this junction tree.
   */
  public String toString() {
    StringBuffer b = new StringBuffer();
    b.append("junction tree over " + varToClusters.keySet() + "\n");
    b.append("with clusters\n");
    ArrayList cl = new ArrayList();
    for (Iterator i = clusters.iterator(); i.hasNext();) {
      Cluster c = (Cluster)i.next();
      cl.add(c);
      b.append("\t" + (cl.size() - 1) + ": " + c.toString() + "\n");
    }
    b.append("and edges (separator sizes)\n");
    for (Iterator i = clusters.iterator(); i.hasNext();) {
      Cluster c = (Cluster)i.next();
      b.append("\t" + cl.indexOf(c) + " ->");
      for (Iterator j = c.edges().iterator(); j.hasNext();) {
	JTEdge e = (JTEdge)j.next();
	Cluster d = (Cluster)e.to;
	b.append(" " + cl.indexOf(d) + "(" + 
		 e.getVariables().size() + ")");
      }
      b.append("\n");
    }
    b.append("and variable index\n");
    for (Iterator i = varToClusters.keySet().iterator(); i.hasNext();) {
      Variable v = (Variable)i.next();
      b.append("\t" + v.toString() + ":");
      for (Iterator j = getClustersWith(v).iterator(); j.hasNext();) 
	b.append(" " + cl.indexOf(j.next()));
      b.append("\n");
    }
    return b.toString();
  }

  /**
   * Tests to see if the junction tree is valid.  A <i>valid</i> junction 
   * tree has three properties: 
   * <ol>
   *   <li><i>potential property</i>: all potentials that have been 
   *       added to this junction tree are covered by at least one 
   *       cluster; (this is not checked)</li>
   *   <li><i>singly-connected property</i>: the graph is a tree; and</li>
   *   <li><i>running intersection property</i>: if <i>C<sub>1</sub></i> and
   *       <i>C<sub>2</sub></i> are clusters that contain the variable
   *       <i>x</i>, then all clusters (and separators) on the (unique) path
   *       between <i>C<sub>1</sub></i> and <i>C<sub>2</sub></i> also 
   *       contain <i>x</i>.</li>
   * </ol>
   * (This method is superfluous, since all public methods of this
   * class preserve the validity of the junction tree.  It is provided
   * as a debugging tool.)
   *
   * @throws InternalError if this junction tree is invalid
   */
  public void checkValid() {
    if (clusters.isEmpty()) return;
    /* Check the singly-connected property by traversing from an
     * arbitrary node and checking to see that each node is visited
     * exactly once.
     */
    Traversal t = new DepthFirstTraversal((Node)clusters.iterator().next());
    while (t.hasNext()) {
      Cluster c = (Cluster)t.next();
      // Check that this cluster is in the index.
      if (!clusters.contains(c))
	  throw new InternalError("Cluster set does not contain " + c
				  + " in " + this);
    }
    if (t.isCyclic())
      throw new InternalError("Junction tree has cycles: " + this);
    if (t.size() < clusters.size()) 
      throw new InternalError("Junction tree is not connected: " + this);
    /* Check the running intersection property.  For each variable and
     * each pair of clusters containing it, find the path between them
     * and check that all clusters and separators on the path contain
     * the variable. */
    for (Iterator i = variables.iterator(); i.hasNext();) {
      Variable var = (Variable)i.next();
      ArrayList clusterList = new ArrayList(getClustersWith(var));
      for (int ci = 0; ci < clusterList.size(); ci++) {
	Cluster c = (Cluster)clusterList.get(ci);
	if (!c.contains(var) || !clusters.contains(c))
	  throw new InternalError("Cluster index invalid: " + this);
	for (int di = ci + 1; di < clusterList.size(); di++) {
	  Cluster d = (Cluster)clusterList.get(di);
	  List path = 
	    c.getShortestPath(new NodeFilter(d) {
		public boolean satisfies(Node n) { return data.equals(n); }});
	  for (Iterator j = path.iterator(); j.hasNext();) {
	    JTEdge e = (JTEdge)j.next();
	    Cluster k = (Cluster)e.to;
	    Cluster h = (Cluster)e.from;
	    /*
	    if (k.getVariables().containsAll(h.getVariables()) ||
		h.getVariables().containsAll(k.getVariables()))
	      throw new InternalError("Nonmaximal clusters: " + this);
	    */
	    if (!k.contains(var))
	      throw new InternalError("Cluster violates " + 
				      "running-intersection property: " +
				      this);
	    if (!clusters.contains(k))
	      throw new InternalError("Cluster list invalid: " + this);
	    if (!e.getPotential().getVariables().contains(var)) 
	      throw new InternalError("Separator on " + e + " violates " + 
				      "running-intersection property for " + 
				      var + " in " + this);
	  }
	}
      }
    }
  }

  /**
   * Tests to see if the junction tree is consistent.  A junction tree
   * is <i>consistent</i> if each separator potential agrees with (the
   * appropriate marginal of) its incident clusters' potentials.
   * (This method is provided as a debugging tool.)
   */
  public boolean consistent() {
    if (clusters.isEmpty()) return true;
    Traversal t = new DepthFirstTraversal((Node)clusters.iterator().next());
    t.next();
    while (t.hasNext()) {
      t.next();
      JTEdge e = (JTEdge)t.edge();
      if (!e.consistent()) return false;
    }
    return true;
  }

  /////////////////////////////////////////////////////////////////////

  /**
   * Returns a directed representation of the junction tree.  An array
   * of integer parent indices into the iteration order of the set
   * returned by {@link JunctionTree#getClusters()}.  The cluster
   * whose parent index is <code>0</code> is the root.
   *
   * <p>This function is provided to facilitate the plotting of
   * junction trees in Matlab.</p>
   *
   * @return An array of integer parent indices into the iteration
   *         order of the set returned by
   *         {@link JunctionTree#getClusters()} 
   */
  public int[] parents() {
    int[] parents = new int[clusters.size()];
    if (clusters.isEmpty()) return parents;
    for (Traversal t = new BreadthFirstTraversal((Node)clusters.get(0));
	 t.hasNext();) {
      Node n = t.next();
      Edge e = t.edge();
      if (e == null)
	parents[clusters.indexOf(n)] = 0;
      else
	parents[clusters.indexOf(n)] = clusters.indexOf(e.from) + 1;
    }
    return parents;
  }

  /**
   * A cluster of a junction tree.
   */
  public class Cluster extends Node {

    /**
     * The potential of this cluster.
     */
    protected Gaussian p;

    /**
     * Creates an empty cluster.
     */
    public Cluster() {
      this(new Gaussian(false));
    }

    /**
     * Creates a cluster with the supplied potential.
     * 
     * @param p a potential
     */
    public Cluster(Gaussian p) {
      this.p = p;
    }

    /**
     * Copy constructor.  This cluster is initialized to have the same
     * potential as the supplied cluster, but to have no edges.
     *
     * @param c a cluster
     */
    public Cluster(Cluster c) {
      this(c.getPotential());
    }

    /**
     * Returns the junction tree this cluster resides in.
     *
     * @return the junction tree this cluster resides in.
     */
    public JunctionTree getJunctionTree() {
      return JunctionTree.this;
    }

    /**
     * Renames a variable in this cluster.
     *
     * @param var   the original variable
     * @param subst the variable substituted for the original variable
     * @throws IllegalArgumentException if <code>var</code> and
     *         <code>subst</code> have differing dimensions
     */
    public void rename(Variable var, Variable subst) {
      p.rename(var, subst);
      for (Iterator i = edges().iterator(); i.hasNext();) {
	JTEdge e = (JTEdge)i.next();
	Cluster c = (Cluster)e.to;
	if (c.contains(var)) e.rename(var, subst);
      }
    }

    /**
     * Returns this cluster's potential.
     */
    public Gaussian getPotential() {
      return p;
    }

    /**
     * Returns the dimension of this cluster, which is the dimension
     * of its potential.
     *
     * @return the dimension of this cluster
     * @see Gaussian#getDimension()
     */
    public int getDimension() {
      return p.getDimension();
    }

    /**
     * Returns the size of this cluster, which is the size of its
     * potential.
     *
     * @return the size of this cluster
     * @see Gaussian#getSize()
     */
    public int getSize() {
      return p.getSize();
    }

    /**
     * Returns an unmodifiable set of this cluster's variables.
     *
     * @return an unmodifiable set of this cluster's variables
     */
    public Set getVariables() {
      return p.getVariables();
    }

    /**
     * Returns <code>true</code> if this cluster is a leaf.
     *
     * @return <code>true</code> if this cluster is a leaf.
     */
    public boolean isLeaf() {
      return (edges.size() == 1);
    }

    /**
     * Returns <code>true</code> if this cluster contains the supplied
     * variable.
     *
     * @return <code>true</code> if this cluster contains <code>v</code>
     */
    public boolean contains(Variable v) {
      return getVariables().contains(v);
    }

    /**
     * Returns <code>true</code> if this cluster contains the supplied
     * set of variables.
     *
     * @param vars a set of {@link Variable Variable}s
     * @return <code>true</code> if this cluster contains the
     *         variables in <code>vars</code>
     */
    public boolean contains(Set vars) {
      return getVariables().containsAll(vars);
    }

    /**
     * Extends this cluster to contain the supplied set of {@link
     * Variable Variable}s.  The potential and the cluster's identifiers
     * are both updated.
     *
     * @param vars a set of {@link Variable Variable}s
     */
    public void extend(Set vars) {
      p.extend(vars);
    }

    /**
     * Marginalizes out a set of variables from this cluster.
     *
     * @param vars a set of {@link Variable Variable}s
     */
    public void marginalizeOut(Set vars) {
      p.marginalizeOut(vars);
    }

    /**
     * Searches for an adjacent cluster that whose variables are all
     * contained in this cluster.  If no such cluster is found,
     * <code>null</code> is returned.
     *
     * @return an adjacent cluster whose variables are all
     *         contained in this cluster, or <code>null</code> if no such 
     *         cluster exists
     */
    public Cluster getSubsumedNeighbor() {
      for (Iterator i = edges().iterator(); i.hasNext();) {
	JTEdge e = (JTEdge)i.next();
	Cluster c = (Cluster)e.to;
	if (contains(c.getVariables())) return c;
      }
      return null;
    }

    /**
     * Searches for an adjacent cluster that contains all the
     * variables contained in this cluster.  If no such cluster is
     * found, <code>null</code> is returned.
     *
     * @return an adjacent cluster that contains all the variables 
     *         in this cluster, or <code>null</code> if no such 
     *         cluster exists
     */
    public Cluster getSubsumingNeighbor() {
      for (Iterator i = edges().iterator(); i.hasNext();) {
	JTEdge e = (JTEdge)i.next();
	Cluster c = (Cluster)e.to;
	if (c.contains(getVariables())) return c;
      }
      return null;
    }

    /**
     * Collects evidence to this cluster.  The neighbors of this
     * cluster (except <code>parent</code>)are asked to collect
     * evidence (with this node as their parent) and this cluster
     * passes a flow to <code>parent</code>.
     *
     * @param local  if <code>true</code>, then the children are not
     *               asked to collect evidence
     * @param parent a neighbor of this cluster that is its parent in
     *               the flow-passing tree, or <code>null</code> if this 
     *               cluster is the root of the flow-passing tree
     */
    public void collectEvidence(Cluster parent, boolean local) {
      JTEdge parentEdge = null;
      for (Iterator i = edges().iterator(); i.hasNext();) {
	JTEdge e = (JTEdge)i.next();
	if (e.to.equals(parent)) continue;
	if (!local)
	  ((Cluster)e.to).collectEvidence(this, false);
	e.getReverse().passFlow();
      }
    }

    /**
     * Distributes evidence from this cluster.  A flow is passed to
     * each neighbor (except <code>parent</code>) and then
     * <code>distributeEvidence</code> is called recursively on all
     * neighbors (except <code>parent</code>).  
     * 
     * @param parent a neighbor of this cluster that is its parent in
     *               the flow-passing tree, or <code>null</code> if this 
     *               cluster is the root of the flow-passing tree
     * @see JunctionTree.Cluster#distributeEvidence(JunctionTree.Cluster,double)
     */
    public void distributeEvidence(Cluster parent) {
      for (Iterator i = edges().iterator(); i.hasNext();) {
	JTEdge e = (JTEdge)i.next();
	if (e.to.equals(parent)) continue;
	e.passFlow();
	((Cluster)e.to).distributeEvidence(this);
      }
    }

    /**
     * Distributes evidence from this cluster so long as the
     * propagated messages are having a significant effect on the
     * belief state.  The significance of a change caused by passing a
     * message is measured in terms of the relative entropy (in nats,
     * i.e., natural logarithmic units) from the old separator
     * potential to its new value after the message is passed.  If a
     * message is passed and the significance is less than a
     * prescribed threshold, then <code>distributeEvidence</code> is
     * not recursively called on the target cluster.
     * 
     * @param parent the parent of this cluster in the flow-passing
     *               tree, or <code>null</code> if this cluster is the root
     * @param thresh the relative entropy (in nats) required from the old
     *               separator potential to the new separator potential
     *               in order to continue distributing evidence
     * @see Gaussian#kl(Gaussian)
     * @see JunctionTree.Cluster#distributeEvidence(JunctionTree.Cluster) 
     */
    public void distributeEvidence(Cluster parent, double thresh) {
      for (Iterator i = edges().iterator(); i.hasNext();) {
	JTEdge e = (JTEdge)i.next();
	if (e.to.equals(parent)) continue;
	Gaussian old = new Gaussian(e.sep);
	e.passFlow();
	if (e.sep.kl(old) >= thresh)
	  ((Cluster)e.to).distributeEvidence(this, thresh);
      }
    }

    /**
     * Returns a string representation of this cluster.
     */
    public String toString() {
      return p.getVariables().toString();
    }
  }

  /**
   * A (directed) edge from one cluster to another in a junction tree.
   * Each (undirected) edge of the junction tree is composed of two
   * directed edges which share the same potential object.
   */
  protected class JTEdge extends Edge {

    /**
     * The separator potential.  This object is shared with the
     * reverse edge associated with this edge, and is thus declared
     * <code>final</code> (so neither edge can alter this sharing).
     */
    private final Gaussian sep;

    /**
     * Constructor.
     *
     * @param from the cluster from which this edge emanates
     * @param to the cluster on which this edge terminates
     * @param p the separator potential of this edge (which is shared
     *          with this edge's complement)
     */
    protected JTEdge(Cluster from, Cluster to, Gaussian p) {
      super(from, to);
      this.sep = p;
    }

    /**
     * Returns the separator potential.
     *
     * @return this edge's separator potential
     */
    public Gaussian getPotential() {
      return sep;
    }

    /**
     * Returns an unmodifiable set of this edge's separator variables.
     *
     * @param an unmodifiable set of this edge's separator variables
     */
    public Set getVariables() {
      return sep.getVariables();
    }


    /**
     * Returns the complement of this edge, i.e., the edge that
     * emanates from this edge's terminus and terminates on this
     * edge's origin (and which shares this edge's separator
     * potential).
     *
     * @return the complement of this edge
     */
    public JTEdge getReverse() {
      return (JTEdge)to.getEdge(from);
    }

    /**
     * Renames a variable in the separator of this edge.
     *
     * @param var   the original variable
     * @param subst the variable substituted for the original variable
     * @throws IllegalArgumentException if <code>var</code> and
     *         <code>subst</code> have differing dimensions
     */
    public void rename(Variable var, Variable subst) {
      sep.rename(var, subst);
    }

    /**
     * Passes a flow along this edge.  The separator potential and the
     * cluster potential of this edge's terminus are updated.
     */
    public void passFlow() {
      if (TJTF.blather) TJTF.say("Passing flow on " + this);
      Gaussian from = ((Cluster)this.from).getPotential();
      Gaussian to = ((Cluster)this.to).getPotential();
      Gaussian fromMarg = from.marginalize(sep.getVariables(), false);
      // This ordering avoids non PSD precision matrices.
      to.times(fromMarg, true);
      to.div(sep, true);
      // It is crucial that we use 'set' here so that this edge's
      // reverse also notices the separator update.
      sep.set(fromMarg);
    }

    /**
     * Returns <code>true</code> iff this edge is consistent.  A
     * junction tree edge is consistent iff the separator marginals of
     * its cluster potentials are equal.
     */
    public boolean consistent() {
      Gaussian from = ((Cluster)this.from).getPotential();
      Gaussian to = ((Cluster)this.to).getPotential();
      // Compute the new separator potential.
      Set vars = sep.getVariables();
      Gaussian fromMarginal = from.marginalize(vars, false);
      Gaussian toMarginal = to.marginalize(vars, false);
      if (!fromMarginal.equals(toMarginal)) {
	if (TJTF.verbose) TJTF.say("Inconsistent edge: " + this);
	return false;
      }
      return true;
    }
  }

  /**
   * Adds a new pair of edges between the supplied clusters.  The
   * separator is initialized to unity.
   *
   * @param c1 a cluster
   * @param c2 another cluster
   */
  protected void link(Cluster c1, Cluster c2) {
    Gaussian sep = new Gaussian(c1.getPotential(), c2.getPotential());
    JTEdge e12 = new JTEdge(c1, c2, sep);
    JTEdge e21 = new JTEdge(c2, c1, sep);
  }
}
