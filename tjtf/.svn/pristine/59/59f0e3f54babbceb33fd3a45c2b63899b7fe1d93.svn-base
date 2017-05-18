package javaslam.slam;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import Jama.*;
import javaslam.util.*;
import javaslam.prob.*;
import javaslam.filter.*;

/**
 *
 * @author Mark A. Paskin
 * @version $Revision: 1.2 $ ($Date: 2003/02/07 20:06:31 $) 
 */
public class LGSLAMFilterCanvas extends Canvas {

  /**
   * The underlying filter.
   */
  protected LGSLAMFilter filter;

  /**
   * Constructor.
   *
   * @param filter  the linear-Gaussian SLAM filter
   */
  public LGSLAMFilterCanvas(LGSLAMFilter filter) {
    this.filter = filter;
  }

  /**
   *
   */
  public void paint(Graphics g) {
    Graphics2D g2 = (Graphics2D)g;
    g2.setStroke(new BasicStroke(2.0f));
    Gaussian p = filter.getRobotMarginal();
    g2.setPaint(Color.red);
    p.paint(g2, 0.95);
    Gaussian[] q = filter.getLandmarkMarginals(null);
    if (q != null) {
      g2.setPaint(Color.green);
      for (int i = 0; i < q.length; i++)
	q[i].paint(g2, 0.95);
    }
  }

  public static java.awt.Frame frame(LGSLAMFilterCanvas c) {
    javax.swing.JFrame f = new javax.swing.JFrame();
    f.addWindowListener(new WindowAdapter() {
	public void windowClosing(WindowEvent e) {System.exit(0);}
      });
    f.getContentPane().add(c);
    f.pack();
    f.setSize(new Dimension(640, 480));
    f.show();
    return f;
  }

  public static LGSLAMFilter f = null;
  
  public static void main(String s[]) {
    double[] mu = new double[2];
    mu[0] = 100;
    mu[1] = 200;
    double[][] sigma = new double[2][];
    sigma[0] = new double[2];
    sigma[1] = new double[2];
    sigma[0][0] = 10.0d;
    sigma[1][1] = 20.0d;
    sigma[1][0] = 8.0d;
    sigma[0][1] = 8.0d;
    f = new KalmanSLAMFilter(mu, sigma);
    frame(new LGSLAMFilterCanvas(f));
  }
}
