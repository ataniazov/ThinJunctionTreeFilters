% SLAMPROB - Description of a SLAM problem structure.
%
% This file describes SLAM problem structures, which are used in a
% number of functions.  For example, SLAM2DPROB creates a SLAM problem
% structure, and SLAM_SIM, SLAM_PROF, SLAM_PLOT, and SLAM_ERR accept a
% SLAM problem structure as an argument.  A SLAM problem structure has
% the following fields:
% 
% dr       - the dimension of the robot state
% dl       - the dimension of each landmark state
% T        - the number of time steps of the problem
% xfun     - a function handle to a function of the form
%                                xn = xfun(x, v, ...)
%            which computes the state of the robot at the next time
%            step as a (possibly non-linear) function of the current
%            robot state x and noise input v.  xfun is assumed to be
%            vectorized so that x and v can be matrices whose columns
%            are individual values (in which case its output must be a
%            matrix whose columns are individual values).
% G        - the positive definite covariance matrix of the state
%            evolution noise variable (v in the description of xfun
%            above) 
% ofun     - a function handle to a function of the form
%                                om = ofun(x, s, ...)
%            which computes the (possibly non-linear) odometry
%            measurement a robot in state x would obtain with
%            white noise s. ofun is assumed to be vectorized so
%            that x and s can be matrices whose columns are
%            individual values (in which case its output must be a
%            matrix whose columns are individual values).
% om       - a matrix with T columns such that odo(:, i) is the
%            odometry measurement at time i
% oC       - the positive definite covariance matrix of the odometry
%            noise variable (s in the description of ofun above)
% lfun     - a function handle to a function of the form
%                            lm = lfun(x, l, w, ...)
%            which computes the (possibly non-linear) measurement a
%            robot in state x would obtain of a landmark in state l
%            with noise w.  lfun is assumed to be vectorized so that x,
%            l, and w can be matrices whose columns are individual
%            values (in which case its output must be a matrix whose
%            columns are individual values).
% ym       - a cell vector such that ym{i} is a matrix whose
%            columns are landmark measurements
% yC       - a positive definite covariance matrix of the landmark
%            measurement noise variable (w in the description of
%            lfun above)
% ilfun    - a function handle to a function of the form
%                              l = ilfun(x, w, lm, ...)
%            which computes the inverse of lfun, i.e., given a robot in
%            state x obtains a measurement lm with noise w, this
%            function computes the state of the landmark l.  ilfun is
%            assumed to be vectorized so that x and v can be matrices
%            whose columns are individual values (in which case its
%            output must be a matrix whose columns are individual
%            values).
%
% The following fields are optional:
%
% path     - a dr-by-T matrix whose columns are the actual states
%            of the robot
% lm       - an dl-by-N matrix whose columns are the states of the
%            landmarks in the problem
% yid      - a cell vector such that yid{i} is a vector of landmark
%            identification numbers where yid{i}(j) is the ID of
%            the landmark that generated measurement ym{i}(:, j)
% largs    - a cell vector such that largs{i} is a cell vector of
%            auxiliary arguments passed to each invocation of lfun
%            at timestep i; this can include a control
%            vector to model active perception, for example
% ilargs   - a cell vector such that ilargs{i} is a cell vector of
%            auxiliary arguments passed to each invocation of ilfun
%            at timestep i; this can include a control
%            vector to model active perception, for example
% xargs    - a cell vector such that xargs{i} is a cell vector of
%            auxiliary arguments passed to the invocation of xfun
%            at timestep i; this can include a control
%            vector, for example
% oargs    - a cell vector such that oargs{i} is a cell vector of
%            auxiliary arguments passed to the invocation of ofun
%            at timestep i; this can include a control
%            vector, for example
%
% The following are examples of functions that generate SLAM problem
% structures: 
%
%   SLAM2DPROB
