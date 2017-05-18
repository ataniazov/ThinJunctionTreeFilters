package javaslam.tjt.graph;

import java.util.*;

/**
 * An iterator that traverses the nodes of a graph.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.3 $ ($Date: 2002/10/14 02:24:11 $)
 */
public abstract class Traversal {

  /**
   * If <code>true</code>, then this iterator has traversed a node
   * more than once.
   */
  protected boolean cycle;
    
  /**
   * A mapping from each node to the {@link Edge Edge} used to arrive
   * at that node (most recently).
   */
  protected HashMap edges;

  /**
   * The most recently returned node, or <code>null</code> if {@link
   * #next() next()} has not yet been called.
   */
  protected Node current;

  /**
   * Constructor.  Implementations of this class should invoke this
   * constructor, and then enqueue the <code>Node</code> parameter.
   *
   * @param n       the node from which to start the traversal; this will be
   *                the first node returned by the iterator
   */
  public Traversal(Node n) {
    this.cycle = false;
    this.current = null;
    edges = new HashMap();
    edges.put(n, null);
  }

  /**
   * Returns <code>true</code> if the iteration has more elements.
   * (In other words, returns true if next would return an element
   * rather than throwing an exception.)
   */
  public abstract boolean hasNext();

  /**
   * Returns <code>true</code> if the traversal has visited a node
   * more than once, implying the existence of a cycle if the graph is
   * undirected.
   */
  public boolean isCyclic() {
    return cycle;
  }

  /**
   * Returns the number of nodes traversed thus far.
   */
  public int size() {
    return edges.keySet().size();
  }

  /**
   * Enqueues a node for later visit.
   */
  protected abstract void enqueue(Node n);

  /**
   * Dequeues a node to visit.
   */
  protected abstract Node dequeue();

  /**
   * Returns the next {@link Node Node} in the iteration.
   */
  public Node next() {
    /* This code is complicated by the fact that the graph can have
     * undirected edges, which can be traversed in both directions.
     * The solution is to record the edge traversed to reach each node
     * in the breadth-first search, and never traverse edges from
     * child back to parent.
     */
    current = dequeue();
    for (Iterator i = current.edges().iterator(); i.hasNext();) {
      Edge e = (Edge)i.next();
      Node c = e.to;
      /* If n and c are connected by an undirected edge, then we may
       * have previously traversed to n from c.  Check to see if this
       * has happened.
       */
      Edge f = (Edge)edges.get(current);
      if ((f != null) && f.from.equals(c)) continue;
      /* If there is an edge recorded for c, then we have previously
       * visited c.
       */
      if (edges.containsKey(c)) {
	cycle = true;
	continue;
      }
      // Enqueue c as a node to be visited, and record the edge used
      // to traverse to it.
      enqueue(c);
      edges.put(c, e);
    }
    return current;
  }

  /**
   * Returns the current node without advancing the iterator.  This
   * method will return <code>null</code> if {@link #next() next()}
   * has not yet been called.
   */
  public Node current() {
    return this.current;
  }

  /**
   * Returns the edge traversed to arrive at the current node.  This
   * method will return <code>null</code> if {@link #next() next()}
   * has not yet been called more than once.
   */
  public Edge edge() {
    if (current == null) return null;
    else return (Edge)edges.get(current);
  }

  /**
   * Returns a list of {@link Edge Edge}s that were traversed in order
   * to reach the Node most recently returned by {@link #next()
   * next()}.
   *
   * @return  a list of the {@link Edge Edge} objects that were
   *          traversed in order to reach the Node most recently 
   *          returned by {@link #next() next()}.
   * @throws IllegalStateException if {@link #next() next()} has not
   *                               yet been called.
   */
  public List path() throws IllegalStateException {
    if (current == null) 
	throw new IllegalStateException("Cannot reconstruct path before" + 
					" next() is called");
    LinkedList list = new LinkedList();
    Edge e = (Edge)edges.get(current);
    while (e != null) {
      list.addFirst(e);
      e = (Edge)edges.get(e.from);
    }
    return list;
  }
}

