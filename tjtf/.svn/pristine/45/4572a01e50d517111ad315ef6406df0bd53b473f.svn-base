% MABSTHRESH - Filters matches based upon a cost matrix using an
%              absolute threshold.
%
% Usage:  a = mabsthresh(a, X, t);
%
% Arguments:   
%         a     - an m x 1 assignment vector.  a(i) is the index of
%                 the feature of the second image that was matched
%                 to feature i of the first image.  If feature i
%                 (of the first image) was not matched to any 
%                 feature in the second image, then a(i) is zero.
%         X     - an m x n affinity matrix
%         t     - a threshold value
%
% Returns:
%         a     - an updated version of a in which each element
%                 that signals a match with cost greater than t
%                 has been zeroed out.
%
% See also MATCH.

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

function a = mabsthresh(a, C, t)

% Remove all matches whose affinity is less than the threshold
b = find(a);
if (~isempty(b))
  v = C(sub2ind(size(C), b, a(b)));
  a(b(find(v > t))) = 0;
end
