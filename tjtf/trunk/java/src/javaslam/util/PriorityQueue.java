package javaslam.util;

import java.util.*;

/**
 * A priority queue.
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.4 $ ($Date: 2002/12/28 21:58:35 $)
 */
public class PriorityQueue {
    
  /**
   * An item in a priority queue.
   */
  protected class QueueElement implements Comparable {

    /**
     * The actual item.
     */
    public Object item;

    /**
     * The priority of the item.
     */
    public Comparable priority;

    /**
     * Constructor.
     */
    public QueueElement(Object item, Comparable priority) {
      this.item = item;
      this.priority = priority;
    }

    /**
     * Compares this object with the specified object for order.
     * Returns a negative integer, zero, or a positive integer as this
     * object is less than, equal to, or greater than the specified
     * object.  The order of two queue elements is exactly the
     * order of their respective priorities.
     *
     * @param o the Object to be compared. 
     * @return a negative integer, zero, or a positive integer as this
     *         object is less than, equal to, or greater than the specified
     *         object.
     * @throws ClassCastException if the specified object's type prevents 
     *                            it from being compared to this Object.
     */
    public int compareTo(Object o) {
      return priority.compareTo(((QueueElement)o).priority);
    }

    /**
     * Returns <code>true</code> if two elements are the same, which
     * is true only when the corresponding items are equal.
     *
     * @param o the Object to be compared. 
     * @return <code>true</code> if two elements are the same.
     */
    public boolean equals(Object o) {
      if (!(o instanceof QueueElement)) return false;
      return item.equals(((QueueElement)o).item);
    }

    /**
     * Returns the hash code of this element's item.
     *
     * @return the hash code of this element's item.
     */
    public int hashCode() {
      return item.hashCode();
    }

    /**
     * Returns a string representation of this queue element.
     */
    public String toString() {
      return item.toString() + "(" + priority.toString() + ")";
    }
  }

  /**
   * An iterator over queue items.  This wraps an iterator over
   * elements and gives access to their corresponding items.
   */
  protected class QueueItemIterator implements Iterator {
    /**
     * The underlying iterator over elements.
     */
    protected Iterator i;

    /**
     * Constructor.
     *
     * @param i an iterator over 
     *          {@link PriorityQueue.QueueElement QueueElement} objects.
     */
    public QueueItemIterator(Iterator i) {
      this.i = i;
    }

    /**
     * Returns <code>true</code> if this iterator has more items.
     */
    public boolean hasNext() {
      return i.hasNext();
    }

    /**
     * Returns the next item in this iteration.
     */
    public Object next() {
      return ((QueueElement)i.next()).item;
    }

    /**
     * This operation is unsupported.
     *
     * @throws UnsupportedOperationException
     */
    public void remove() {
      throw new UnsupportedOperationException();
    }
  }
  
  /**
   * A list of {@link QueueElement QueueElement} objects ordered by
   * their priorities.
   */
  protected ArrayList elements;

  /**
   * Maps items to their elements.
   */
  protected Map itemsToElements;

  /**
   * Constructor.
   */
  public PriorityQueue() {
    elements = new ArrayList();
    itemsToElements = new HashMap();
  }
  
  /**
   * Returns <code>true</code> iff this queue is empty.
   */
  public boolean isEmpty() {
    return elements.isEmpty();
  }

  /**
   * Returns an iterator over the items in this priority queue in
   * order of decreasing priority.
   */
  public Iterator iterator() {
    return new QueueItemIterator(elements.iterator());
  }

  /**
   * Returns the number of items in this queue.
   */
  public int size() {
    return elements.size();
  }

  /**
   * Removes all items from this queue.
   */
  public void clear() {
    elements.clear();
    itemsToElements.clear();
  }

  /**
   * Enqueues a new item with the supplied priority.  If the item was
   * already present in the queue, then this method has no effect; in
   * particular, it does not change the priority of <code>item</code>
   * to <code>priority</code>.
   *
   * @param item     the new item to enqueue
   * @param priority the priority of the new item
   * @throws ClassCastException if the supplied priority is not 
   *                            comparable to the priorities of the 
   *                            previously enqueued items
   */
  public void enqueue(Object item, Comparable priority) {
    if (itemsToElements.containsKey(item)) return;
    QueueElement element = new QueueElement(item, priority);
    int idx = Collections.binarySearch(elements, element);
    if (idx < 0) idx = -idx - 1;
    itemsToElements.put(item, element);
    elements.add(idx, element);
  }
   
  /**
   * Dequeues the item with the highest priority.  If this queue is
   * empty, then <code>null</code> is returned.
   *
   * @return the item with the highest priority or <code>null</code>
   *          if this queue is empty
   */
  public Object dequeue() {
    if (elements.isEmpty()) return null;
    QueueElement last = (QueueElement)elements.get(elements.size() - 1);
    if (itemsToElements.remove(last.item) == null)
      throw new Error("Internal inconsistency");
    if (!elements.remove(last))
      throw new Error("Internal inconsistency");
    return last.item;
  }

  /**
   * Removes a particular item from the queue, if it exists.
   *
   * @param item the item to be removed
   */
  public void remove(Object item) {
    QueueElement element = (QueueElement)itemsToElements.get(item);
    if (element == null) return;
    if (itemsToElements.remove(item) == null)
      throw new Error("Internal inconsistency");
    if (!elements.remove(element))
      throw new Error("Internal inconsistency");
  }

  /**
   * Returns the item with the highest priority without removing it
   * from the queue.  If this queue is empty, then <code>null</code>
   * is returned.
   *
   * @return the item with the highest priority or <code>null</code>
   *          if this queue is empty
   */
  public Object peek() {
    if (elements.isEmpty()) return null;
    QueueElement last = (QueueElement)elements.get(elements.size() - 1);
    return last.item;
  }

  /**
   * Returns the priority of the highest priority item in the queue.
   *
   * @return the highest priority of any item in the queue or
   *          <code>null</code> if this queue is empty
   */
  public Comparable highest() {
    if (elements.isEmpty()) return null;
    QueueElement last = (QueueElement)elements.get(elements.size() - 1);
    return last.priority;
  }

  /**
   * Updates the priority of an item already in the queue.  If the
   * item was not in the queue, then it is enqueued with the supplied
   * priority.
   *
   * @param item     the item whose priority is to be changed
   * @param priority the new priority of the item
   * @throws ClassCastException if the supplied priority is not 
   *                            comparable to the priorities of the 
   *                            previously enqueued items
   */
  public void reprioritize(Object item, Comparable priority) {
    remove(item);
    enqueue(item, priority);
  }

  /**
   * Returns a string representation of this priority queue.
   */
  public String toString() {
    return elements.toString();
  }
}
