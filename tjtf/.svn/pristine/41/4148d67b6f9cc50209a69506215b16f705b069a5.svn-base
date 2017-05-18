package javaslam.tjt.graph;

/**
 * A predicate over graph nodes.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.2 $ ($Date: 2002/10/11 16:57:14 $)
 */
public abstract class NodeFilter {

  /**
   * An arbitrary piece of data used by the predicate.
   */
  protected Object data;

  /**
   * Constructor.
   *
   * @param data the local data used by this predicate
   */
  public NodeFilter(Object data) {
    this.data = data;
  }

  /**
   * Returns <code>true</code> if the supplied node satisfies this
   * predicate.
   *
   * @param n the node to which the predicate is applied
   * @return <code>true</code> iff <code>n</code> satisfies this predicate
   */
  public abstract boolean satisfies(Node n);
}
