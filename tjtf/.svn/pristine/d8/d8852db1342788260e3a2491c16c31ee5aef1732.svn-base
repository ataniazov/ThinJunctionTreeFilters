% LGSLAM_MAP - Extracts the map estimate of a linear-Gaussian
%              SLAM filter.
%
% Usage:
% 
%   [x, lm, id] = lgslam_err(f)
%
% Inputs:
%
%   f     - a javaslam.slam.LGSLAMFilter object
%
% Outputs:
%
%   x   - the robot state estimate
%   lm  - a matrix whose columns are landmark state estimates
%   id  - a vector whose entries give the landmark ids
%         corresponding to lm

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
function [x, lm, id] = lgslam_map(bs)
  % Get the marginals over all the variables.
  map = bs.getMarginals([]);
  xp = map.get(bs.getRobotVariable);
  x = xp.getMu([]).getArray;
  % Extract the landmark estimates.
  id = bs.getLandmarkIds;
  lm = zeros(2, length(id));
  for i=1:length(id)
    lp = map.get(bs.getLandmarkVariable(id(i)));
    lm(:, i) = lp.getMu([]).getArray;
  end
