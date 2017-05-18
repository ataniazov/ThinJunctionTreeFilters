package javaslam.tjt.graph;

import java.util.*;

/**
 * An iterator that traverses the nodes of a graph in breadth-first
 * (siblings before children) order.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.2 $ ($Date: 2002/10/11 16:57:14 $)
 */
public class BreadthFirstTraversal extends Traversal {

  /**
   * A queue of nodes to be visited in this search.
   */
  protected LinkedList queue;

  /**
   * Constructor.
   *
   * @param n       the node from which to start the traversal; this will be
   *                the first node returned by the iterator
   */
  public BreadthFirstTraversal(Node n) {
    super(n);
    this.queue = new LinkedList();
    enqueue(n);
  }

  /**
   * Returns <code>true</code> if the iteration has more elements.
   * (In other words, returns true if next would return an element
   * rather than throwing an exception.)
   */
  public boolean hasNext() {
    return !queue.isEmpty();
  }

 /**
  * Enqueues a node for later visit.
  */
  protected void enqueue(Node n) {
    queue.add(n);
  }

  /**
   * Dequeues a node to visit.
   */
  protected Node dequeue() {
      return (Node)queue.removeFirst();
  }
}

