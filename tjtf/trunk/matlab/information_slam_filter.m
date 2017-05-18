% INFORMATION_SLAM_FILTER - Generates a new SLAM filter structure whose
%                           belief state is maintained using the Information
%                           filter.
%
% Usage:
%
%   f = information_slam_filter(m, C)
%
% Inputs: 
%
%   m - the starting state of the robot
%   C - a positive definite matrix giving the covariance
%       (uncertainty) in the robot's initial state
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
function sf = information_slam_filter(m, C)

sf.label = 'Information filter';
sf.bs = javaslam.slam.InformationSLAMFilter(m, C);
sf.ffun = @lgslam;
sf.dafun = @lgslam_da;
sf.pfun = @lgslam_plot;
sf.mfun = @lgslam_map;

return;
