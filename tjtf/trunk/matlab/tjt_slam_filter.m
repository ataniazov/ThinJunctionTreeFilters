% TJT_SLAM_FILTER - Generates a new SLAM filter structure whose
%                   belief state is a thin junction tree.
%
% Usage:
%
%   f = tjt_slam_filter(m, C, width, overlap[, sig])
%
% Inputs: 
%
%   m       - the starting state of the robot
%   C       - a positive definite matrix giving the covariance
%             (uncertainty) in the robot's initial state
%   width   - the maximum cluster size of the junction tree 
%   overlap - the maximum separator size of the junction tree
%   sig     - if non-empty, this is the significance threshold
%             used in adaptive message passing (default: [])
%
% See also: SLAMFILTER

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
function sf = tjt_slam_filter(m, C, width, overlap, significance)


sf.bs = javaslam.slam.TJTSLAMFilter(width, overlap, m, C);
if (nargin > 4)
  sf.label = sprintf('TJTF (k = %d, h = %d, s = %0.2g)', ...
		     width, overlap, significance);
  sf.bs.getJunctionTree.setSignificance(significance);
else
  sf.label = sprintf('TJTF (k = %d, h = %d)', ...
		     width, overlap);
end
sf.ffun = @lgslam;
sf.dafun = @lgslam_da;
sf.pfun = @jt_slam_plot;
sf.mfun = @lgslam_map;

return;
