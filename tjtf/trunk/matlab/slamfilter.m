% SLAMFILTER - Describes a SLAM filter structure.
%
% A SLAM filter structure has the following fields:
% 
%   label: a string giving a short description of the filter
%   bs:    the belief state of the filter represented as a Java
%          object 
%   ffun:  a handle to a function that performs a filtering
%          operation to update the belief state; for a detailed
%          description of its parameters, see LGSLAM or FASTSLAM,
%          which are examples.
%   dafun: a handle to a function that performs data association;
%          for a detailed description of its parameters, see LGSLAM_DA,
%          which is an example.
%   pfun:  a handle to a function that plots the belief state; for
%          a detailed description of its parameters, see
%          KALMAN_SLAM_PLOT, TJT_SLAM_PLOT or FASTSLAM_PLOT, which
%          are examples.
%   mfun:  a handle to a function that computes the estimated
%          map; for a detailed description of its parameters, see 
%          LGSLAM_MAP or FASTSLAM_MAP, which are examples.
%
% The following are examples of functions that generate SLAM filter
% objects: 
%
%   KALMAN_SLAM_FILTER
%   INFORMATION_SLAM_FILTER
%   JT_SLAM_FILTER
%   TJT_SLAM_FILTER
%   FASTSLAM_FILTER

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
