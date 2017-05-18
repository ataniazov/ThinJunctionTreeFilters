% KALMAN_SLAM_PLOT - Plots a Gaussian distribution that is the belief state
%                    of a Kalman SLAM filter.
%
% Inputs:
%
%   f - a javaslam.slam.KalmanSLAMFilter object
%
% Options:
% 
%   colors      - a cell array of color specifiers; colors{i} is the
%                 color in which the marginal of landmark with
%                 identifier i should be plotted; if this is [],
%                 then all landmarks are plotted in black (default: [])
%   robot-color - the color in which the robot's marginal should be
%                 plotted (default: 'm')
%   conf        - the size of the confidence ellipsoids plotted
%                 (default: 0.95)

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

function [] = lgslam_plot(f, varargin)

[colors, ...
 conf, ...
 robot_color, ...
 R, ...
 t, ...
 bw] = process_options(varargin, ...
		       'colors', [], ...
		       'conf', 0.9, ...
		       'robot-color', 'm', ...
		       'rotation', [], ...
		       'translation', [], ...
		       'b&w', 0);

% Extract the marginals of all landmarks.
ids = f.getLandmarkIds;
lm = f.getLandmarkMarginals(ids);

hold on;

% Plot the covariance matrices of all landmarks.
for i=1:length(lm)
  m = lm(i).getMu([]).getArray;
  C = lm(i).getSigma([], []).getArray;
  if (~isempty(R))
    m = R * m(1:2);
    C = R * C(1:2, 1:2) * R';
  end
  if (~isempty(t))
    m = m(1:2) + t;
  end
  % Plot the covariance ellipsoid.
  if (bw)
    color = 'k';
  elseif (~isempty(colors))
    color = colors{ids(i)};
  else
    color = 'g';
  end
  line_width = 0.5;
  plotcov2(m(1:2), C(1:2, 1:2), 'conf', conf, ...
	   'plot-opts', {'Color', color, 'LineWidth', line_width}, ...
	   'num-pts', 20);
end

% Extract the marginals of the robot state.
robot = f.getRobotMarginal;

hold on;

% Plot the covariance matrices of all landmarks.
m = robot.getMu([]).getArray;
C = robot.getSigma([], []).getArray;
if (~isempty(R))
  m = R * m(1:2);
  C = R * C(1:2, 1:2) * R';
end
if (~isempty(t))
  m = m(1:2) + t;
end
% Plot the covariance ellipsoid.
if (bw)
  color = 'k';
else
  color = robot_color;
end
line_width = 0.5;
plotcov2(m(1:2), C(1:2, 1:2), 'conf', conf, ...
	 'plot-opts', {'Color', color, 'LineWidth', line_width}, ...
	 'num-pts', 20);
