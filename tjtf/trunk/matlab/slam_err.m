% SLAM_ERR
%
% Computes error measures for the belief state of a SLAM filter
% operating on a SLAM problem.
%
% Usage:
%
%   slam_err(p, f, t[, OPTIONS])
%
% Inputs:
%
%   p  - a SLAM problem structure (see SLAMPROB)
%   f  - a SLAM filter structure
%   t  - the current time step
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
function [xm_error, lm_avg_error, lm_std_error, lm_max_error] = ...
    slam_err(p, f, t)

% Extract the filter's map and compute the error in robot position.
[x, lm, id] = feval(f.mfun, f.bs);
xm_error = norm(x(1:2) - p.path(1:2, t));
% Align the filter's map to the true map using a rigid transformation.
[c, R, t] = ralign(lm, p.lm(:, id));
lm = R * lm + repmat(t, [1, length(id)]);
% Compute the landmark and some statistics.
lm_err = sqrt(sum((lm - p.lm(:, id)).^2, 1));
lm_avg_error = mean(lm_err);
lm_std_error = std(lm_err);
lm_max_error = max(lm_err);
