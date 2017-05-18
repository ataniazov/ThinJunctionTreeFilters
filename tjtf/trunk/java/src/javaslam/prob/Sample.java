package javaslam.prob;

import java.util.Random;
import javaslam.util.*;

/**
 * A utility class containing methods for sampling from various
 * probability distributions.  The implementations used here are
 * borrowed directly from Tom Minka's Lightspeed library for Matlab.
 *
 * <p>This class records counts of all floating point operations using
 * {@link Flops#count(long)}.</p>
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.4 $ ($Date: 2003/02/07 02:26:06 $) 
 */
public class Sample {

  /**
   * Analogous to <tt>eps</tt> in Matlab.
   */
  public final static double DBL_EPSILON = 2.2204e-16;

  /**
   * Samples from the Gamma(a) distribution.  This implementation is
   * from BUGS.
   */
  public static double gamma(double a, Random r) {
    double c1 = 0.0d, c2 = 0.0d, c3 = 0.0d, c4 = 0.0d, c5 = 0.0d, c6 = 0.0d;
    if (Double.isNaN(a)) return a;
    if (a < DBL_EPSILON) return 0.0d;
    /* If a == 1, then gamma is exponential. */
    /* It is important that random() is never zero. */
    if (Math.abs(a - 1) < DBL_EPSILON) {
      Flops.count(Flops.log() + Flops.rand());
      return -Math.log(r.nextDouble());
    }
    if (a < 1.0d) {
      c6 = Math.exp(Math.log(r.nextDouble()) / a);
      a = a + 1.0d;
      Flops.count(Flops.exp() + Flops.log() + Flops.rand() + 2);
    } else c6 = 1.0d;
    c1 = a - 1.0d;
    c2 = (a - 1.0d / (6.0d * a)) / c1;
    c3 = 2.0d / c1;
    c4 = c3 + 2.0d;
    c5 = 1.0d / Math.sqrt(a);
    c6 *= c1;
    Flops.count(9 + Flops.sqrt());
    while (true) {
      double u1 = 0.0d, u2 = 0.0d, w = 0.0d;
      /* loop until u1 is valid */
      while (true) {
	u1 = r.nextDouble();
	u2 = r.nextDouble();
	Flops.count(2 * Flops.rand());
	if (a > 2.5d) {
	  u1 = u2 + c5 * (1 - 1.86d * u1);
	  Flops.count(4);
	}
	if (!((u1 <= 0) || (u1 >= 1))) break;
      }
      w = c2 * u2 / u1;
      Flops.count(2);
      if (c3 * u1 + w + 1.0d / w <= c4) {
	Flops.count(4);
	return c6 * w;
      }
      if (c3 * Math.log(u1) - Math.log(w) + w < 1.0d) {
	Flops.count(3 + 2 * Flops.log());
	return c6 * w;
      }
    }
  }

  /**
   * Returns a sample from Beta(a, b).
   */
  public static double beta(double a, double b, Random r) {
    double g = gamma(a, r);
    Flops.count(2);
    return g / (g + gamma(b, r));
  }

  /**
   * Returns a sample from Binomial(n, p).
   */
  public static int binomial(int n, double p, Random r) {
    int t = 0;
    if (Double.isNaN(p)) return 0;
    if (p < DBL_EPSILON) return 0;
    if (p >= 1 - DBL_EPSILON) return n;
    if (n < 15) {
      /* Coin flip method. This takes O(n) time. */
      for (int i = 0; i < n; i++) 
	if (r.nextDouble() < p) t++;
      Flops.count(n * Flops.rand());
      return t;
    }
    Flops.count(1);
    if ((double)n * p < 10.0d) {
      /* Waiting time method.  This takes O(np) time. */
      double q = -Math.log(1-p), e = -Math.log(r.nextDouble()), s = 0.0d;
      Flops.count(2 * Flops.log() + Flops.rand());
      t = n;
      for(s = e/(double)t; s <= q; s += e/(double)t) {
	t--;
	if (t == 0) break;
	e = -Math.log(r.nextDouble());
      }
      t = n - t;
      Flops.count(1 + t * (3 + Flops.log() + Flops.rand()));
      return t;
    }
    /* Recursive method.  This makes O(log(log(n))) recursive calls. */
    int i = (int)Math.round(Math.floor(p * (n + 1)));
    Flops.count(1);
    double b = Sample.beta(i, n + 1 - i, r);
    if (b <= p) {
      t = i + Sample.binomial(n - i, (p - b) / (1 - b), r);
      Flops.count(3);
    } else {
      t = i - Sample.binomial(i - 1, (b - p) / b, r);
      Flops.count(2);
    }
    return t;
  }
}
