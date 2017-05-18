% LGSLAM - Performs a filter update for a linear-Gaussian SLAM
%          filter (using the unscented transformation to
%          linearize nonlinear motion and measurement models).
%  
% This function first incorporates all landmark observations, then
% evolves the state according to the robot's dynamics, and then
% performs an odometry update.
%
% Usage:
% 
%   lgslam(f, lfun, ilfun, lm, lC, id, xfun, G, ofun, om, oC[, OPTIONS]);
%
% Inputs:
%
%   f     - a javaslam.slam.LGSLAMFilter object
%   lfun  - a function handle to a function of the form
%                            lm = lfun(x, l, w, ...)
%           which computes the (possibly non-linear) measurement a
%           robot in state x would obtain of a landmark in state l
%           with noise w.  lfun is assumed to be vectorized so that x,
%           l, and w can be matrices whose columns are individual
%           values (in which case its output must be a matrix whose
%           columns are individual values).  If this function is not
%           supplied, the robot is assumed to receive no landmark
%           measurements.
%   ilfun - a function handle to a function of the form
%                              l = ilfun(x, v, lm, ...)
%           which computes the inverse of lfun, i.e., given a robot in
%           state x obtains a measurement lm with noise v, this
%           function computes the state of the landmark l.  ilfun is
%           assumed to be vectorized so that x and v can be matrices
%           whose columns are individual values (in which case its
%           output must be a matrix whose columns are individual
%           values).  If this function is not supplied, new landmarks
%           are not added to the map.
%   lm    - an h x m matrix whose columns are landmark observations
%   lC    - a positive definite (possibly empty) matrix giving the
%           covariance of the white noise supplied to the landmark
%           observation function lfun
%   id    - a data association vector; id(i) gives the landmark
%           identifier of the landmark that generated measurement
%           lm(:, i)
%   xfun  - a function handle to a function of the form
%                                xn = xfun(x, v, ...)
%           which computes the state of the robot at the next time
%           step as a (possibly non-linear) function of the current
%           robot state x and noise input v.  xfun is assumed to be
%           vectorized so that x and v can be matrices whose columns
%           are individual values (in which case its output must be a
%           matrix whose columns are individual values).  If this
%           function is not supplied, the robot state is assumed
%           not to have changed.
%   G     - a positive semi-definite (possibly empty) matrix giving
%           the covariance of the white noise supplied to the state
%           transition function xfun
%   ofun  - a function handle to a function of the form
%                                om = ofun(x, s, ...)
%           which computes the (possibly non-linear) odometry
%           measurement a robot in state x would obtain with
%           white noise s. ofun is assumed to be vectorized so
%           that x and s can be matrices whose columns are
%           individual values (in which case its output must be a
%           matrix whose columns are individual values). If this
%           function is not supplied, the robot is assumed to
%           receive no odometry measurements. 
%   om    - a k x 1 vector giving the odometry measurement; if this
%           is [], then no odometry measurement update is performed
%   oC    - a positive definite (possibly empty) matrix giving
%           the covariance of the white noise supplied to the odometry
%           observation function ofun
%
% Options:
%
%   'largs'   - a cell array of auxiliary arguments passed to each
%               invocation of lfun; this can include a control
%               vector to model active perception, for example; the
%               default is {}.
%   'ilargs'  - a cell array of auxiliary arguments passed to each
%               invocation of ilfun; this can include a control
%               vector to model active perception, for example; the
%               default is {}.
%   'xargs'   - a cell array of auxiliary arguments passed to each
%               invocation of xfun; this can include a control
%               vector, for example; the default is {}.
%   'oargs'   - a cell array of auxiliary arguments passed to each
%               invocation of ofun; this can include a control
%               vector, for example; the default is {}.
%   'verbose' - if nonzero, then this function displays progress
%               messages (default: 0)
%
% This function counts flops.  (See FLOPS.)

% Copyright (C) 2002 Mark A. Paskin
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
% USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lgslam(f, lfun, ilfun, lm, lC, id, xfun, G, ofun, om, oC, ...
		  varargin)

import javaslam.util.*;
import java.util.*;

[xargs, ...
 largs, ...
 ilargs, ...
 oargs] = process_options(varargin, 'xargs', {}, ...
			            'largs', {}, ...
			            'ilargs', {}, ...
			            'oargs', {});

n = f.getNumLandmarks; % The number of landmarks
m = size(lm, 2);       % The number of landmark measurements
dr = f.getRobotVariable.dim;

% Cached marginal of the robot state.
xm = [];
xC = [];

% If there are landmark measurements, integrate them.
if (~isempty(lfun) & (m > 0))
  % Find the landmarks that are currently in the filter.
  known = [];
  unknown = [];
  for k=1:m
    if (f.contains(id(k)))
      known = [known; k];
    else
      unknown = [unknown; k];
    end
  end
  % Integrate the known landmarks.
  if (~isempty(known))
    ps = f.getRobotLandmarkMarginals(id(known));
    % Cache the robot's marginal.
    xm = ps(1).getMu(ListSet(f.x)).getArray;
    xC = ps(1).getSigma(ListSet(f.x), []).getArray;
    for i=1:length(known)
      k = known(i);
      y = lm(:, k);                     % Observation vector
      % This landmark has been observed before.  Get the marginal over
      % the robot and landmark state.
      v = f.getLandmarkVariable(id(k));
      zm = ps(i).getMu(ListSet([f.x, v])).getArray;
      zC = ps(i).getSigma(ListSet([f.x, v]), []).getArray;
      % Perform the update using the Unscented Kalman Filter.
      [a, B, Q] = ut(zm, zC, lC, @lm_obs_model, 'args', {dr, lfun, largs{:}});
      f.measurement(id(k), full(a), full(B(:, 1:f.x.dim)), ...
		    full(B(:, (f.x.dim + 1):end)), full(Q), full(y));
    end
  end
  % Integrate the unknown landmarks
  if (~isempty(unknown))
    if (isempty(xm))
      p = f.getRobotMarginal;
      if p.isCanonical
	p.reparameterize(1);
      end
      % Cache the robot's marginal.
      xm = p.getMu(ListSet(f.x)).getArray;
      xC = p.getSigma(ListSet(f.x), []).getArray;
    end
    for i=1:length(unknown)
      k = unknown(i);
      y = lm(:, k);                     % Observation vector
      % This is a first-time landmark observation.  In order to get a
      % good linearization of the observation function, we have to first
      % get a rough estimate of the landmark's (as yet unobserved)
      % state.  To do this, we use (a linearization of) the observation
      % function's inverse.
      [a, B, Q, ym, yC, CC] = ut(xm, xC, lC, ilfun, ...
				 'args', {y, ilargs{:}});    
      % Now ym and yC give an estimate of the landmark's state.  Use
      % it to linearize the forwards measurement model.
      zm = [xm; ym];
      zC = [xC, CC; CC', yC];
      % Perform the update using the Unscented Kalman Filter.
      [a, B, Q] = ut(zm, zC, lC, @lm_obs_model, 'args', {dr, lfun, largs{:}});
      f.measurement(id(k), full(a), full(B(:, 1:f.x.dim)), ...
		    full(B(:, (f.x.dim + 1):end)), full(Q), full(y));
    end    
  end
end

% Time update.
if (~isempty(xfun))
  % If it has not yet been obtained, compute and cache the robot's
  % state marginal.  Note that if the value is cached, the cache
  % may not be accurate: incorporating the landmark observations
  % above could have significantly altered the marginal.  However,
  % since we're using it for linearization only, it is (often) more
  % efficient to ignore this discrepancy.
  if (isempty(xm))
    % Get the marginal over the robot's position.
    p = f.getRobotMarginal;
    if p.isCanonical
      p.reparameterize(1);
    end
    xm = p.getMu([]).getArray;
    xC = p.getSigma([], []).getArray;
  end
  % Linearize the motion model.
  [a, B, Q] = ut(xm, xC, G, xfun, 'args', xargs);
  % Apply the update to the filter.
  f.motion(full(a), full(B), full(Q));
  % Update the cached marginal of the robot state efficiently.
  xm = a + B * xm;
  xC = B * xC * B' + Q;
  d = length(a);
  Flops.count(d + Flops.mult(d, d, 1) + 2 * Flops.mult(d, d, d) + d^2);
end

% Perform the odometry update
if (~isempty(om) & ~isempty(ofun))
  if (isempty(xm))
    % Get the new marginal over the robot's position.
    p = f.getRobotMarginal;
    if p.isCanonical
      p.reparameterize(1);
    end
    xm = p.getMu([]).getArray;
    xC = p.getSigma([], []).getArray;
  end
  % Linearize the odometry model.
  [a, B, Q] = ut(xm, xC, oC, ofun, 'args', oargs);
  % Apply the update to the filter.
  f.odometry(full(a), full(B), full(Q), full(om));
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Applies the landmark observation function to the sub-vectors
% that represent the robot's and landmark's states.
function obs = lm_obs_model(state, noise, dr, lfun, varargin)
  obs = feval(lfun, state(1:dr, :), state((dr+1):end, :), noise, varargin{:});
  return;
