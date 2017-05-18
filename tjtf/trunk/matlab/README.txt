============================================================================
SETUP.
============================================================================

These Matlab files are a convenient interface for the Java library
containing the implementation of thin junction tree filters (TJTF).
They make use of the fact that Matlab has an internal Java virtual
machine (JVM).  To learn more about this, see 

  http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_external/ch_java2.shtml

To have the TJTF java code loaded into the Matlab JVM, your Matlab
classpath.txt file must include references to the root of the javaslam
package, and also to Jama-1.0.1.jar (both of which are contained in
the bin subdirectory of tjtf-java.zip).

============================================================================
A SIMPLE EXAMPLE.

Once your classpath.txt file is set up and the Matlab files in this
distribution are in your Matlab path, you should be able to execute
this simple example.  Once you've run the example, you can trace
through the code of SLAM_SIM to see how it all works.
============================================================================

% First, create a simple two-dimensional SLAM problem.
p = slam2dprob('num-landmarks', 100, 'side', 20, 'shape', 'square');

% Let's take a look at what we just made.  
slam2dplot(p);

% The landmarks are shown as black dots.  The robot's path is the blue
% line.  The red line is the integrated odometry signal, and the green
% line is the integrated control signal.  The cyan dots are the
% landmark measurements (relative to the robot's true position).

% Now let's build a Kalman filter for this problem.  Start by making
% an initial covariance for the robot state, which in this case has five
% dimensions.
xC = 1e-5 * eye(5); 

% Now create a filter with the initial covariance and the true
% starting state of the robot.
f = kalman_slam_filter(p.path(:, 1), xC);

% Okay, let's watch the Kalman filter do its thing!
slam_sim(p, f);

============================================================================
FILES.

Below is a categorization of the files, along with brief descriptions.
Most of the files are very well documented and commented, so it should
be easy to figure out what they do and how they fit together.  (Some
of the comments and documentation are slightly out of date.)
============================================================================

Building a SLAM problem:
-----------------------

slamprob.m			# describes a SLAM problem structure
  slam2dprob.m			# generates an example SLAM problem
    h2rot.m			# converts angle to rotation matrix
    slam2dplot.m		# plots an example SLAM problem

Building a SLAM filter:
----------------------

slamfilter.m			# describes a SLAM filter structure
  kalman_slam_filter.m		# builds a Kalman SLAM filter 
  information_slam_filter.m	# builds an information SLAM filter
  jt_slam_filter.m		# builds a junction tree SLAM filter
    tjt_slam_filter.m		# builds a thin junction tree SLAM filter


Things to do with a SLAM filter:
-------------------------------

slam_sim.m			# visualizes a SLAM filter on a SLAM problem
  slam_plot.m			# visualizes the belief state of a filter
slam_prof.m			# compares several SLAM filters on a problem
  slam_err.m			# computes the estimation error of a filter

Internal SLAM code:
----------------------

lgslam.m			# performs a linear-Gaussian filter update 
lgslam_da.m			# performs maximum likelihood data association
  match.m			# wrapper for the hungarian method
  hungarian.m			# weighted bipartite matching code
  mahal2conf.m			# used for gating in data association
  mabsthresh.m			# used for gating in data association
lgslam_map.m			# extracts map from a linear-Gaussian filter
lgslam_plot.m			# plots the belief state of a LG filter
jt_slam_plot.m			# plots the beliefs of a junction tree filter


Plotting code:
-------------

conf2mahal.m			# computes confidence ellipse size
plotcov2.m			# plots a 2D confidence ellipse
plotcov3.m			# plots a 3D confidence bubble
checkpsd.m			# checks for positive-definite covariance


Miscellaneous code:
------------------

process_options.m		# argument/option processing used by all files
ralign.m			# Rigid alignment of two 2D data sets










