package javaslam.tjt.graph;

import java.util.*;

/**
 * A node in a graph.  The {@link Node#equals(Object) equals} and
 * {@link Node#hashCode hashCode} methods are implemented exactly as
 * in the {@link Object Object} class.  Thus, it is not possible to
 * have two representations of the same node.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.3 $ ($Date: 2002/10/18 23:33:56 $)
 */
public class Node extends Object {

  /**
   * A mapping from neighbor {@link Node}s to the {@link Edge}s to them. 
   */
  protected Map edges;

  /**
   * Default constructor.
   */
  public Node() {
    this.edges = new HashMap();
  }

  /**
   * Adds an edge from this node to another node.  If the edge's
   * {@link Edge#from from} field is not <code>this</code>, then it is
   * set to <code>this</code>.  
   *
   * @param e the edge to add to this node's list of outward edges
   * @throws IllegalArgumentException if the edge's
   *         {@link Edge#to to} field is equal to <code>this</code>.
   */
  public void addEdge(Edge e) {
    if (this.equals(e.to))
      throw new IllegalArgumentException("Cannot add self-loops");
    e.from = this;
    edges.put(e.to, e);
  }

  /**
   * Changes an edge emanating from this node so that it emanates from
   * the supplied node.  If the terminus of the edge is the same as
   * the supplied node, then the edge is simply deleted.  This method
   * has no effect if the supplied edge does not emanate from this
   * node.
   *
   * @param e the edge to swing
   * @param n the new terminus of <code>e</code>
   */
  public void swingEdge(Edge e, Node n) {
    if (!e.from.equals(this)) return;
    this.removeEdge(e);
    if (!e.to.equals(n))
      n.addEdge(e);
    // If the edge is undirected, we must redirect its reverse.
    Edge f = e.to.getEdge(this);
    if (f == null) return;
    // Remove the edge, redirect it, and add it back.
    e.to.removeEdge(f);
    f.to = n;
    if (!e.to.equals(n))
      e.to.addEdge(f);
  }

  /**
   * Changes the origin of all (directed and undirected) edges
   * emanating from this node to the supplied node.  Edges between
   * this node and the supplied node are dropped.
   *
   * @param n the new origin of all edges emanating from this node
   */
  public void swingEdgesTo(Node n) {
    while (!edges.isEmpty()) {
      Edge e = (Edge)edges().iterator().next();
      swingEdge(e, n);
    }
  }

  /**
   * Gets the edge from this node to the supplied node; if no such edge
   * exists, <code>null</code> is returned.
   *
   * @param n the neighbor of this node
   * @return the edge from this node to <code>n</code>, or
   *         <code>null</code> if <code>n</code> is not a neighbor of this
   *         node
   */
  public Edge getEdge(Node n) {
    return (Edge)edges.get(n);
  }

  /**
   * Removes the supplied edge from this node to the supplied node (if
   * such an edge exists).
   *
   * @param e an edge emanating from this node
   */
  public void removeEdge(Edge e) {
    edges.remove(e.to);
  }

  /**
   * Returns (an unmodifiable view of) the set of {@link Edge Edge}s
   * incident from this node.
   *
   * @return (an unmodifiable view of) the set of {@link Edge Edge}s
   *         emanating from this node
   */
  public Collection edges() {
    return Collections.unmodifiableCollection(edges.values());
  }

  /**
   * Computes the shortest path from this node to another node
   * satisfying a particular property.  The path is returned as a
   * <code>List</code> of {@link Edge Edge} objects.  If there are no
   * nodes reachable from this node that satisfy the filter, then
   * <code>null</code> is returned.
   *
   * @param f a filter identifying nodes that have a particular property
   * @return  a path to the node closest to this one (in path length) 
   *          that satisfies the filter (or <code>null</code> if no node
   *          satisfying the property is reachable from this node).
   */
  public List getShortestPath(NodeFilter f) {
    for (Traversal t = new BreadthFirstTraversal(this); t.hasNext();) {
      Node n = t.next();
      if (f.satisfies(n)) 
	try { return t.path(); } catch (IllegalStateException X) {
	  throw new Error("Could not get path.");
	}
    }
    return null;
  }
    
  /**
   * Finalizes the definition of <code>equals</code> so that two nodes
   * are equal iff they are the same object.
   *
   * @param o another object
   * @return <code>true</code> iff <code>o == this</code>
   */
  public final boolean equals(Object o) {
    return super.equals(o);
  }

  /**
   * Finalizes the definition of <code>hashCode</code> to be the same as
   * that of {@link Object#hashCode hashCode}.
   */
  public final int hashCode() {
    return super.hashCode();
  }
}

