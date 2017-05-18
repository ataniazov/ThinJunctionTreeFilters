% LGSLAM_DA - Maximum-likelihood data association for
%             linear-Gaussian SLAM filters with gating. 
%
% We build an m x n cost matrix L such that L(i, j) is the smallest
% confidence interval containing measurement i (under the state
% distribution of landmark j).  We then use L as the input to an
% matching problem; i.e., we search for the assignment of measurements
% to landmarks of minimum cost.  Then matches whose smallest
% confidence interval exceeds the supplied threshold are rematched to
% new landmarks.
%
% Usage:
% 
%   id = lgslam_da(f, lfun, lm, lC[, OPTIONS]);
%
% Inputs:
%
%   f    - a javaslam.slam.LGSLAMFilter object
%   lm   - an h x m matrix whose columns are landmark observations
%   lC   - a positive definite (or empty) 
%          matrix giving the covariance of the white noise supplied to 
%          the landmark observation function lfun
%   lfun - a function handle to a function of the form
%                           lm = lfun(x, l, w, ...)
%          which computes the (possibly non-linear) measurement a
%          robot in state x would obtain of a landmark in state
%          l with noise w.
%
% Options:
%
%   'args'  - a cell array of auxiliary arguments passed to each
%             invocation of lfun; this can include a control
%             vector to model active perception, for example; the
%             default is {}.
%   'conf'  - a probability threshold giving the minimum likelihood that
%             is permissible for a measurement-landmark assignment; the
%             measurement must lie within the landmark's observation
%             confidence ellipse with the supplied probability
%             (default: 0.95)
%
% Outputs:
%
%   id   - an m x 1 vector of integers identifying which landmark
%          each measurement was assigned to; the i'th measurement
%          was assigned to landmark id(i)
%
% This function counts flops.  (See FLOPS.)
%

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

function id = lgslam_da(f, lfun, lm, lC, varargin)

import javaslam.slam.*;
import javaslam.util.*;
import java.util.*;

[args, conf] = process_options(varargin, 'args', {}, 'conf', 0.95);

n = f.getNumLandmarks;  % The number of current landmarks
m = size(lm, 2);       % The number of landmark measurements
x = f.getRobotVariable; % The variable representing the robot state

if ((n > 0) & (m > 0))
  % If there are assignments and measurements, then do data
  % association.
  L = zeros(m, n);
  % In order to populate the matching matrix L, we must consider each
  % landmark/observation pair separately (because the landmark's
  % observation density depends upon its state estimate as well as the
  % observation noise, which varies across measurements).  This can be
  % very slow; a faster implementation would use an index to limit
  % consideration to landmarks sufficiently close (in Euclidean
  % distance) to the measurements.

  % Get the robot-landmark joint marginals efficiently.
  ids = f.getLandmarkIds;
  ps = f.getRobotLandmarkMarginals(ids);
  for i=1:n
    vs = ListSet([x f.getLandmarkVariable(ids(i))]);
    um(:, i) = ps(i).getMu(vs).getArray;
    uC(:, :, i) = ps(i).getSigma(vs, []).getArray;
  end
    
  % Populate the matching matrix.
  for lmid=1:n    
    for mid=1:m
      % Use the Unscented Kalman Filter to get a predictive density
      % P(z(t + 1) | z(1:t), lm(t + 1) = lmid).
      y = lm(:, mid);                    % Observation vector
      R = lC;                            % Observation noise
      [a, B, G, ym, yC] = ut(um(:, lmid), uC(:, :, lmid), ...
			     R, @lm_obs_model, ...
			     'args', {x.dim, lfun, args{:}});
      % Using this density, compute the confidence interval defined by
      % the measurement.
      inn = y - ym;
      tmp = length(y);
      L(mid, lmid) = mahal2conf(inn' * (yC \ inn), tmp);
      Flops.count(tmp + Flops.solve(tmp, tmp, 1, 1) + ...
		  Flops.mult(1, tmp, 1));
    end
  end
  % Remove all landmarks that don't have a measurement associated
  % with the threshold probability.
  z = find(min(L, [], 1) <= conf);
  if (isempty(z))
    id = (n + 1):(n + m);
  else
    L = L(:, z);
    % Match the measurements to landmarks
    id = match(L);
    % Threshold the matches
    id = mabsthresh(id, L, conf);
    % Map back to landmark indexes.
    id(find(id)) = z(id(find(id)));
    % Map from landmark indexes to landmark IDs.
    id(find(id)) = ids(id(find(id)));
    % Give fresh indexes to landmarks that were not assigned.
    k = find(~id);
    id(k) = (n + 1):(n + length(k));
  end
elseif ((n == 0) & (m > 0))
  % If there are measurements but no landmarks, then each is
  % allocated a new landmark.
  id = (1:m)';
else
  % There are no measurements; we can return early.
  id = zeros(1, 0);
  return;
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Applies the landmark observation function to the sub-vectors
% that represent the robot's and landmark's states.
function obs = lm_obs_model(state, noise, dr, lfun, varargin)
  obs = feval(lfun, state(1:dr, :), state((dr+1):end, :), noise, varargin{:});
  return;
