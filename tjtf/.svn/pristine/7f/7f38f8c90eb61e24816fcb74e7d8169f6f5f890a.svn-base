% SLAM2DPROB - Generates a 2D SLAM problem structure.  
%
% A world consisting of a square plane with point landmarks is
% constructed; the robot travels along a prescribed path (subject to
% velocity limits), observing the landmarks that lie within a
% cone-shaped visual field.
%
% Options:
%
%   'seed'                 - the seed used to initialize both RAND and
%                            RANDN so that the simulation generated is
%                            a deterministic function of this
%                            functions inputs; if this is [], then the
%                            random number generator is not seeded
%   'shape'                - the shape of the robot path; this can be
%                            'square', 'circle', 'line', or
%                            'switchbacks' (default: 'square')
%   'side'                 - the length of the side of the square region
%                            to be mapped (in meters)
%   'num-landmarks'        - the number of landmarks (default: 1000)
%   'random-landmarks'     - a flag; if this is nonzero, the
%                            landmark locations are sampled
%                            uniformly on the square region; if
%                            not, they are spaced on a grid
%                            (default: 1)
%   'max-range'            - defines the length of the visual cone
%                            (in meters; default: 10 m)
%   'max-bearing'          - defines half of the arc of the visual cone
%                            (in radians; default: pi / 2)
%   'num-obs'              - if nonzero, then the maximum range is
%                            chosen dynamically to ensure that this
%                            num-obs landmark observations are made
%                            at each time step (when this is
%                            possible; default: [])
%   'max-trans-ctrl'       - the maximum translation velocity of
%                            the robot (default: 0.5 meters/second)
%   'max-rot-ctrl'         - the maximum rotational velocity of the
%                            robot (default: 30*pi/180 radians/second)
%   'ctrl-rel-trans-noise' - the variance of the relative noise in
%                            the robot's translation control signal
%                            (default: 0.03^2)
%   'ctrl-rel-rot-noise'   - the variance of the relative noise in
%                            the robot's rotation control signal
%                            (default: 0.05^2)
%   'ctrl-abs-trans-noise' - the variance of the absolute noise in
%                            the robot's translation control signal
%                            (default: (0.02 meters/sec)^2)
%   'ctrl-abs-rot-noise'   - the variance of the absolute noise in
%                            the robot's rotation control signal
%                            (default: (pi/180 radians/sec)^2)
%   'odo-trans-noise'      - the variance of the absolute noise in
%                            the robot's translation odometry
%                            (default: (0.05 meters/sec)^2)
%   'odo-rot-noise'        - the variance of the absolute noise in
%                            the robot's rotation odometry
%                            (default: (pi/180 radians/sec)^2)
%   'bearing-noise'        - the variance of the absolute noise in
%                            the robot's bearing measurements
%                            (default: (2 * pi/180 radians/sec)^2)
%   'rel-range-noise'      - the variance of the relative noise in
%                            the robot's range measurements
%                            (default: (0.1)^2)
%   'abs-range-noise'      - the variance of the absolute noise in
%                            the robot's range measurements
%                            (default: (0.5 meters)^2)
%
% See also: SLAMPROB, SLAM2DPLOT
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

function p = slam2dprob(varargin)

[seed, ...
 shape, ...
 side, ...
 random_lm, ...
 num_obs, ...
 max_range, ...
 max_bearing, ...
 max_trans_ctrl, ...
 max_rot_ctrl, ...
 ctrl_rel_trans_noise, ...
 ctrl_rel_rot_noise, ...
 ctrl_abs_trans_noise, ...
 ctrl_abs_rot_noise, ...
 odo_trans_noise, ...
 odo_rot_noise, ...
 bearing_noise, ...
 rel_range_noise, ...
 abs_range_noise, ...
 n] = process_options(varargin, 'seed', [], ...
			        'shape', 'square', ...
			        'side', 100, ...
			        'random-landmarks', 1, ...
		                'num-obs', [], ...
			        'max-range', 10, ...
			        'max-bearing', pi / 2, ...
			        'max-trans-ctrl', 0.5, ...
			        'max-rot-ctrl', (30 * pi/180), ...
			        'ctrl-rel-trans-noise', 0.03^2, ...
			        'ctrl-rel-rot-noise', 0.05^2, ...
			        'ctrl-abs-trans-noise', 0.02^2, ...
			        'ctrl-abs-rot-noise', (1 * pi / 180)^2, ...
			        'odo-trans-noise', (0.05)^2, ...
			        'odo-rot-noise', (1 * pi / 180)^2, ...
			        'bearing-noise', (2 * pi / 180)^2, ...
			        'rel-range-noise', 0.1^2, ...
			        'abs-range-noise', 0.5^2, ...
			        'num-landmarks', 1000);

if (~isempty(seed))
  % It is necessary to seed these two random number generators
  % independently.
  rand('state', seed);
  randn('state', seed);
end

% The dimensions of the robot and landmark states
p.dr = 5;  % X = [x, y, h, t, r]
p.dl = 2;  % L = [x, y]

% Robot state evolution noise covariance.
p.G = diag([ctrl_abs_trans_noise, ctrl_rel_trans_noise, ...
	    ctrl_abs_rot_noise, ctrl_rel_rot_noise]);

% Odometric noise covariance.
p.oC = diag([odo_trans_noise, odo_rot_noise]);

% Landmark measurement noise covariance.
p.yC = diag([bearing_noise, rel_range_noise, abs_range_noise]); 

% Generate the desired pose of the robot at every time step.
if (strcmpi(shape, 'square'))
  % The nubmer of steps necessary to traverse a side and make a
  % 90-degree turn.
  side_steps = ceil(side / max_trans_ctrl);
  rot_steps = ceil((pi / 2) / max_rot_ctrl);
  % First side.
  traj = [linspace(0, side, side_steps); zeros(2, side_steps)];
  % First rotation.
  traj = [traj, [repmat(traj(1:2, end), [1 rot_steps]); ...
		 linspace(0, -pi / 2, rot_steps)]];
  % Second side.
  traj = [traj, [repmat(traj(1, end), [1 side_steps]); ...
		 linspace(0, side, side_steps); ...
		 repmat(traj(3, end), [1 side_steps])]];
  % Second rotation.
  traj = [traj, [repmat(traj(1:2, end), [1 rot_steps]); ...
		 linspace(-pi / 2, -pi, rot_steps)]];
  % Third side.
  traj = [traj, [linspace(side, 0, side_steps); ...
		 repmat(traj(2, end), [1 side_steps]); ...
		 repmat(traj(3, end), [1 side_steps])]];
  % Third rotation.
  traj = [traj, [repmat(traj(1:2, end), [1 rot_steps]); ...
		 linspace(-pi, -3 * pi / 2, rot_steps)]];
  % Fourth side.
  traj = [traj, [repmat(traj(1, end), [1 side_steps]); ...
		 linspace(side, 0, side_steps); ...
		 repmat(traj(3, end), [1 side_steps])]];
elseif (strcmpi(shape, 'circle'))
  steps = ceil(side * pi / max_trans_ctrl);
  rot_v = 2 * pi / steps;
  trans_v = side * pi / steps;
  traj = zeros(3, steps);
  traj(3, :) = linspace(0, 2 * pi, steps);
  for i=2:steps
    traj(1:2, i) = traj(1:2, i - 1) + h2rot(traj(3, i)) * [trans_v; 0];
  end
elseif (strcmpi(shape, 'line'))
  steps = ceil(side / max_trans_ctrl);
  traj = zeros(3, steps);
  traj(1, :) = linspace(0, side, steps);
elseif strcmpi(shape, 'switchbacks')
  % The nubmer of steps necessary to traverse a side and make a
  % 90-degree turn.
  side_steps = ceil(side / max_trans_ctrl);
  rot_steps = ceil((pi / 2) / max_rot_ctrl);
  % Choose the number of switchbacks so that the field of view
  % overlaps on consecutive passes.
  sb_side = max_range * sin(max_bearing);
  num_sb = round(side / (sb_side * 2));
  sb_steps = ceil(sb_side / max_trans_ctrl);
  % Compute the trajectory.
  traj = [0; 0; 0];
  for i=1:num_sb
    % First long side.
    traj = [traj, [linspace(0, side, side_steps); ...
		   repmat(traj(2, end), [1 side_steps]); ...
		   repmat(traj(3, end), [1 side_steps])]];
    % First rotation.
    traj = [traj, [repmat(traj(1:2, end), [1 rot_steps]); ...
		   linspace(0, -pi / 2, rot_steps)]];
    % First small side.
    traj = [traj, [repmat(traj(1, end), [1 sb_steps]); ...
		   traj(2, end) + linspace(0, sb_side, sb_steps); ...
		   repmat(traj(3, end), [1 sb_steps])]];
    % Second rotation.
    traj = [traj, [repmat(traj(1:2, end), [1 rot_steps]); ...
		   linspace(-pi / 2, -pi, rot_steps)]];
    % Second long side.
    traj = [traj, [linspace(side, 0, side_steps); ...
		   repmat(traj(2, end), [1 side_steps]); ...
		   repmat(traj(3, end), [1 side_steps])]];
    % Third rotation.
    traj = [traj, [repmat(traj(1:2, end), [1 rot_steps]); ...
		   linspace(-pi, -pi / 2, rot_steps)]];
    % Second small side.
    traj = [traj, [repmat(traj(1, end), [1 sb_steps]); ...
		   traj(2, end) + linspace(0, sb_side, sb_steps); ...
		   repmat(traj(3, end), [1 sb_steps])]];
    % Final rotation.
    traj = [traj, [repmat(traj(1:2, end), [1 rot_steps]); ...
		   linspace(-pi / 2, 0, rot_steps)]];
  end
  % One more long side.
  traj = [traj, [linspace(0, side, side_steps); ...
		 repmat(traj(2, end), [1 side_steps]); ...
		 repmat(traj(3, end), [1 side_steps])]];
else
  error(sprintf('Unknown shape: %s', shape));
end
p.T = size(traj, 2);

% Create a robot state trajectory from this pose trajectory.
p.path = [traj; zeros(2, p.T)];
% Reverse-engineer the translational velocity.
p.path(4, 1:(p.T - 1)) = ...
    sqrt((traj(1, 2:p.T) - traj(1, 1:(p.T - 1))).^2 + ...
	 (traj(2, 2:p.T) - traj(2, 1:(p.T - 1))).^2);
p.path(4, p.T) = 0;
% Reverse-engineer the rotational velocity.
p.path(5, 1:(p.T - 1)) = ...
    mod(traj(3, 2:p.T) - traj(3, 1:(p.T - 1)) + pi, 2 * pi) - pi;
p.path(5, p.T) = 0;

% Now 'invert' the control noise to obtain an omniscient controller
% that makes the robot traverse the desired path, in spite of
% control noise.
ctrl = zeros(2, p.T);
% Compute the desired change in velocity.
del = p.path(4:5, 2:p.T);
% Compute the noise at every time step.
noise = mvnrnd(zeros(size(p.G, 1), 1), p.G, p.T - 1)';
% Subtract out the additive noise.
del = del - noise([1 3], :);
% Divide out the multiplicative noise.
del = del ./ (1 + noise([2 4], :));
ctrl(:, 1:(p.T - 1)) = del;
if 1
  % Check to make sure we've inverted the dynamics correctly.
  for i=1:(p.T - 1)
    p.path(:, i + 1) = rse(p.path(:, i), noise(:, i), ctrl(:, i));
  end
end

% Set up the arguments to the system functions.
p.largs = cell(1, p.T); [p.largs{:}] = deal({});
p.oargs = cell(1, p.T); [p.oargs{:}] = deal({});
p.ilargs = cell(1, p.T); [p.ilargs{:}] = deal({});
p.xargs = cell(1, p.T);
for k=1:p.T
  p.xargs{k} = { ctrl(:, k) };
end

% Compute the (noisy) odometry from the path.
p.om = odo_obs(p.path, mvnrnd([0; 0], p.oC, p.T)');

% Compute a region with twice the are of the bounding box of the path.
top = max(p.path(1:2, :), [], 2) + repmat(max_range, [2 1]);
bottom = min(p.path(1:2, :), [], 2) - repmat(max_range, [2 1]);

p.ym = cell(1, p.T);
p.yid = cell(1, p.T);
p.lm = [];

if (n > 0)
  if (random_lm)
    % Populate the world with landmarks scattered uniformly over the
    % plane.
    p.lm = (rand(2, n) .* repmat(top - bottom, [1, n])) + ...
	   repmat(bottom, [1, n]);
  else
    % Place the landmarks in a uniform grid.
    n = floor(sqrt(n));
    [lmx, lmy] = ...
	meshgrid(linspace(bottom(1), top(1), n), ...
		 linspace(bottom(2), top(2), n));
    p.lm = [lmx(:)'; lmy(:)'];
    n = n^2;
  end
  % Compute the (noisy) landmark measurements.
  for k=1:p.T
    % Compute the coordinates of the landmarks in the robot's current pose.
    rview = lm_obs(repmat(p.path(:, k), [1 n]), p.lm, zeros(3, n));
    
    % Compute the robot-centric range and bearing of the landmarks.
    [b, r] = cart2pol(rview(1, :), rview(2, :));
    
    % Make those in a small, forward cone visible.
    if (isempty(num_obs))
      p.yid{k} = find((abs(b) <= max_bearing) & (r <= max_range));
    else
      cidx = find(abs(b) <= max_bearing);
      if (length(cidx) <= num_obs)
	p.yid{k} = cidx;
      else
	[tmp, sidx] = sort(r(cidx));
	p.yid{k} = cidx(sidx(1:num_obs));
      end
    end
    nobs = length(p.yid{k});
    
    % Generate the observations.
    p.ym{k} = lm_obs(repmat(p.path(:, k), [1 nobs]), ...
		     p.lm(:, p.yid{k}), ...
		     mvnrnd(zeros(size(p.yC, 1), 1), p.yC, nobs)');
  end 
end

p.xfun = @rse;
p.ofun = @odo_obs;
p.lfun = @lm_obs;
p.ilfun = @lm_obs_inv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function is the robot's state evolution function.  The control
% signal gives a rotation and displacement, and the noise vector is
% added into the control signal.  First the robot is displaced, and
% then it is rotated.
function ns = rse(s, noise, control)
  [m, n] = size(s);
  ns = zeros(m, n);
  % Set the new position and heading to the previous position and
  % heading plus the previous velocities.
  tr = zeros(2, n);  tr(1, :) = s(4, :);
  ns(1, :) = s(1, :) + s(4, :) .* cos(s(3, :));
  ns(2, :) = s(2, :) + s(4, :) .* sin(s(3, :));
  ns(3, :) = s(3, :) + s(5, :);
  % The position and rotation derivatives are supplied by the control
  % (plus absolute and relative noise).
  ns(4:5, :) = repmat(control, [1 n]) .* (1 + noise([2 4], :)) ...
      + noise([1 3], :);
  return;

% This function is the robot's odometry measurement.
function obs = odo_obs(s, noise)
  % The odometry consists of the derivatives of the translation and
  % rotation.
  obs = s(4:5, :) + noise;
  return;

% This function is the robot's landmark measurement model.  The
% landmark measurements are given in Cartesian coordinates
% (relative to the robot), but the noise terms are expressed in
% polar (bearing, range) coordinates.
function obs = lm_obs(rs, ls, noise)
  obs = zeros(size(ls));
  % Compute the landmark's location in the robot's (Cartesian)
  % coordinate frame.
  y = ls(1:2, :) - rs(1:2, :);
  % Compute the observation in polar coordinates.
  obs(1, :) = atan2(y(2, :), y(1, :)) - rs(3, :);  % bearing
  obs(2, :) = sqrt(y(1, :).^2 + y(2, :).^2);       % range
  % Inject the position noise.
  obs(1, :) = obs(1, :) + noise(1, :);                      
  obs(2, :) = obs(2, :) .* (noise(2, :) + 1) + noise(3, :);
  % Convert to Cartesian coordinates.
  [obs(1, :), obs(2, :)] = pol2cart(obs(1, :), obs(2, :));
  % Inject the appearance noise.
  obs(3:end, :) = ls(3:end, :) + noise(4:end, :);
  return;

% This function is the inverse of the robot's landmark measurement
% model.  The landmark measurements are given in Cartesian coordinates
% (relative to the robot).
function ls = lm_obs_inv(rs, noise, obs)
  % Compute the number of points.
  n = size(rs, 2);
  % Allocate space for the landmark states.
  ls = zeros(2, n);
  % Replicate the observation to be the same size as rs and noise.
  if (size(obs, 2) ~= n)
    obs = repmat(obs, [1 n]);
  end
  % Convert the observation to polar coordinates.
  [obs(1, :), obs(2, :)] = cart2pol(obs(1, :), obs(2, :));
  % Remove the observation noise.
  obs(1, :) = obs(1, :) - noise(1, :);                      % bearing
  obs(2, :) = obs(2, :) ./ (noise(2, :) + 1) - noise(3, :); % range
  % Compute the landmark's location in the global (Cartesian)
  % coordinate frame.
  ls(1, :) = rs(1, :) + obs(2, :) .* cos(obs(1, :) + rs(3, :));
  ls(2, :) = rs(2, :) + obs(2, :) .* sin(obs(1, :) + rs(3, :));
  return;
  