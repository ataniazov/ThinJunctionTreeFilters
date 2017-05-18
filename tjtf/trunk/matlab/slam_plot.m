% SLAM_PLOT
%
% Visualizes the belief state of a SLAM filter operating on a SLAM
% problem.
%
% Usage:
%
%   slam_plot(p, f[, OPTIONS])
%
% Inputs:
%
%   p  - a SLAM problem structure (see SLAMPROB)
%   f  - a SLAM filter structure
% 
% Options:
%
%   'plot-opts'  - the options passed to the plotting method
%   'resolution' - if non-empty, this must be a two vector 
%                  [height width] giving the resolution
%                  of the figure window (default: [])
%   'time'       - if this is non-empty, the current location of
%                  the robot is plotted 
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

function [] = slam_plot(p, f, varargin)

[plot_opts, ...
 resolution, ...
 time, ...
 align] = process_options(varargin, ...
			  'plot-opts', {}, ...
			  'resolution', [], ...
			  'time', [], ...
			  'align', 0);

% Clear the figure and bring it forward.
figure(gcf);
clf;
set(gcf, 'Units', 'pixels');
set(gcf, 'DoubleBuffer', 'on');
if (~isempty(resolution))
  pos = get(gcf, 'Position');
  set(gcf, 'Position', [pos(1), pos(2), resolution(1), resolution(2)]);
end
bottom = min([p.lm, p.path(1:2, :)], [], 2);
top = max([p.lm, p.path(1:2, :)], [], 2);
spread = (top - bottom);
top = top + spread * 0.25;
bottom = bottom - spread * 0.25;
set(gca, 'XLim', [bottom(1) top(1)]);
set(gca, 'YLim', [bottom(2) top(2)]);
set(gca, 'DataAspectRatio', [1 1 1]);
set(gca, 'Visible', 'off');
set(gcf, 'Color', [1 1 1]);
hold on;
    
if (~isempty(p.lm))
  % Plot the true landmarks' locations.
  plot(p.lm(1, :), p.lm(2, :), 'k.');
end
    
if (~isempty(time))
  % Plot the robot's true location.
  plot(p.path(1, time), p.path(2, time), 'ro', ...
       'MarkerSize', 6, 'LineWidth', 2);
end

% If an alignment is requested, compute the best rigid transformation.
if (align)
  [x, lm, id] = feval(f.mfun, f.bs);
  [c, R, t] = ralign(lm, p.lm(:, id));
  plot_opts = {plot_opts{:}, 'rotation', R, 'translation', t};
end

% Plot the belief state.
feval(f.pfun, f.bs, plot_opts{:});
    
drawnow;

