package javaslam.tjt.graph;

/**
 * A directed edge in a graph.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.2 $ ($Date: 2002/10/11 16:57:14 $)
 */
public class Edge {

  /**
   * The start of this edge.
   */
  public Node from;

  /**
   * The end of this edge.
   */
  public Node to;

  /**
   * Default constructor.  The edge is built and added to the
   * <code>from</code> node's list of incident edges.
   */
  public Edge(Node from, Node to) {
    this.from = from;
    this.to = to;
    from.addEdge(this);
  }

  /**
   * Redefines <code>equals</code> so that two edges are equal if
   * their termini are the equal.
   */
  public boolean equals(Object o) {
    return ((o instanceof Edge) && 
	    (((Edge)o).to.equals(to)) &&
	    (((Edge)o).from.equals(from)));
  }

  /**
   * Redefines <code>hashCode</code> so the hash code of an edge is
   * the sum of the hash codes of its nodes.
   */
  public int hashCode() {
    return to.hashCode() + from.hashCode();
  }

  /**
   * Returns a string representation of this edge.
   */
  public String toString() {
    return from.toString() + " -> " + to.toString();
  }
}
