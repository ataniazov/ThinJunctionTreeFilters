% H2ROT - Converts a heading angles to a 2 by 2 rotation matrix.
%
% Usage:
%   A = pry2rot(h)
% 
% Inputs:
%   h - the heading angle (clockwise about the z-axis from the
%       positive x-axis in radians)
%
% Output:
%   A - a 2 by 2 matrix such that a point w = [x y]', when subjected to
%       the rotation identified by h is at A * w.
%
% This function is vectorized so that h may be a column vector.  In
% this case, A is a 3D matrix whose pages (slices) are rotation
% matrices corresponding to the elements of h.

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
function A = h2rot(h);

A = zeros(2, 2, length(h));
A(1, 1, :) = cos(h);
A(1, 2, :) = sin(h);
A(2, 1, :) = -sin(h);
A(2, 2, :) = cos(h);


