% JT_SLAM_PLOT - Plots a junction tree that is the belief state
%                of a SLAM filter.
%
% Inputs:
%
%   f - a javaslam.slam.JTSLAMFilter object
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
%   type        - the method of displaying dependencies:
%                   edges: edges are plotted between each pair of
%                          variables in the same cluster; this is the
%                          Gaussian graphical model
%                   jtree:  the junction tree is plotted, with clusters
%                           represented by square nodes,
%                           cluster-variable edges representing
%                           membership, and cluster-cluster edges
%                           representing adjacency in the junction tree
%                   submap: for each cluster, a the smallest box that
%                           encloses the cluster's landmarks is plotted;
%                           the box is blue if the cluster also contains
%                           the robot state variable
%                   none:   dependencies are not displayed
%                 (The default is 'jtree'.)
%   bold-edge-color,
%   edge-color  - color specifiers; if 'type' is 'edges' or
%                 'jtree', then regular edges are plotted using
%                 edge-color (default: [0.5 0.5 0.5]) and special
%                 edges (involving the robot) are plotted using
%                 bold-edge-color (default: [0.5 0.5 0.5])

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

function [] = jt_slam_plot(f, varargin)

[colors, ...
 conf, ...
 type, ...
 robot_color, ...
 edge_color, ...
 bold_edge_color, ...
 R, ...
 t] = process_options(varargin, ...
		      'colors', [], ...
		      'conf', 0.9, ...
		      'type', 'jtree', ...
		      'robot-color', 'm', ...
		      'edge-color', [0.5 0.5 0.5], ...
		      'bold-edge-color', [0.8 0.8 0.8], ...
		      'rotation', [], ...
		      'translation', []);

n = f.getNumLandmarks;

% Extract the means and covariances of all state variables.
map = f.getMarginals([]);

hold on;
if strcmp(type, 'edges')
  % Plot the edges of the graphical model.
  clusters = f.getJunctionTree.getClusters;
  c_iter = clusters.iterator;
  while (c_iter.hasNext)
    c = c_iter.next;
    v_iter1 = c.getVariables.iterator;
    while (v_iter1.hasNext)
      v1 = v_iter1.next;
      v_iter2 = c.getVariables.iterator;
      while (v_iter2.hasNext)
	v2 = v_iter2.next;
	if (~v1.equals(v2))
	  a = map.get(v1).getMu([]).getArray;
	  b = map.get(v2).getMu([]).getArray;
	  if (~isempty(R))
	    a = R * a(1:2);
	    b = R * b(1:2);
	  end
	  if (~isempty(t))
	    a = a(1:2) + t;
	    b = b(1:2) + t;
	  end
	  % Color edges from the robot state dark gray, and all
	  % inter-landmark edges light gray.
	  if (v1.equals(f.x) | v2.equals(f.x))
	    color = bold_edge_color;
	  else
	    color = edge_color;
	  end
	  plot([a(1), b(1)], [a(2), b(2)], 'Color', color);
	end
      end
    end
  end
elseif strcmp(type, 'jtree')
  cluster_set = f.getJunctionTree.getClusters;
  parents = double(f.getJunctionTree.parents);
  clusters = cell(cluster_set.size, 1); % Cluster 
  c_iter = cluster_set.iterator;
  k = 0;
  % Plot lines from the cluster nodes to the variables.
  while (c_iter.hasNext)
    k = k + 1;
    clusters{k} = c_iter.next;
    % Compute the location of the cluster node.
    cloc{k} = zeros(2, 1);
    v_iter = clusters{k}.getVariables.iterator;
    while (v_iter.hasNext)
      v = v_iter.next;
      a = map.get(v).getMu([]).getArray;
      if (~isempty(R))
	a = R * a(1:2);
      end
      if (~isempty(t))
	a = a(1:2) + t;
      end
      cloc{k} = cloc{k} + a(1:2);
    end
    cloc{k} = cloc{k} / clusters{k}.getVariables.size;
    v_iter = clusters{k}.getVariables.iterator;
    while (v_iter.hasNext)
      v = v_iter.next;
      a = map.get(v).getMu([]).getArray;
      if (~isempty(R))
	a = R * a(1:2);
      end
      if (~isempty(t))
	a = a(1:2) + t;
      end
      plot([a(1), cloc{k}(1)], [a(2), cloc{k}(2)], 'Color', edge_color);
    end
  end
  % Plot the edge between the cluster nodes.
  c_iter = cluster_set.iterator;
  for k=1:cluster_set.size
    a = cloc{k};
    if (parents(k))
      b = cloc{parents(k)};
      plot([a(1), b(1)], [a(2), b(2)], '-', 'Color', bold_edge_color);
    end
  end
  % Plot the cluster nodes.
  for k=1:cluster_set.size
    a = cloc{k};
    plot(a(1), a(2), 's', ...
	 'MarkerSize', 6, 'MarkerFaceColor', 'k', ...
	 'MarkerEdgeColor', 'k');
  end
elseif strcmp(type, 'submap')
  % Plot the submaps.
  cluster_set = f.getJunctionTree.getClusters;
  clusters = cell(cluster_set.size, 1); % Cluster 
  cmax = cell(cluster_set.size, 1);     % Location of first corner 
  cmin = cell(cluster_set.size, 1);     % Location of second corner 
  c_iter = cluster_set.iterator;
  k = 0;
  while (c_iter.hasNext)
    k = k + 1;
    clusters{k} = c_iter.next;
    % Compute the corners of each map.
    cmax{k} = -Inf * ones(2, 1);
    cmin{k} = Inf * ones(2, 1);
    v_iter = clusters{k}.getVariables.iterator;
    while (v_iter.hasNext)
      v = v_iter.next;
      if (v.equals(f.x)) 
	continue;
      end
      a = map.get(v).getMu([]).getArray;
      if (~isempty(R))
	a = R * a(1:2);
      end
      if (~isempty(t))
	a = a(1:2) + t;
      end
      u = find(a(1:2) < cmin{k});
      cmin{k}(u) = a(u);
      u = find(a(1:2) > cmax{k});
      cmax{k}(u) = a(u);
    end
  end
  % Plot the maps.
  for k=1:cluster_set.size
    if (clusters{k}.contains(f.x))
      color = robot_color;
    else
      color = 'k';
    end
    plot([cmin{k}(1); cmax{k}(1)], ...
	 [cmin{k}(2); cmin{k}(2)], 'Color', color);
    plot([cmin{k}(1); cmin{k}(1)], ...
	 [cmin{k}(2); cmax{k}(2)], 'Color', color);
    plot([cmin{k}(1); cmin{k}(1)], ...
	 [cmin{k}(2); cmin{k}(2)], 'Color', color);
    plot([cmax{k}(1); cmax{k}(1)], ...
	 [cmin{k}(2); cmax{k}(2)], 'Color', color);
    plot([cmax{k}(1); cmax{k}(1)], ...
	 [cmin{k}(2); cmin{k}(2)], 'Color', color);
    plot([cmin{k}(1); cmax{k}(1)], ...
	 [cmax{k}(2); cmax{k}(2)], 'Color', color);
    plot([cmin{k}(1); cmin{k}(1)], ...
	 [cmax{k}(2); cmax{k}(2)], 'Color', color);
    plot([cmin{k}(1); cmax{k}(1)], ...
	 [cmin{k}(2); cmin{k}(2)], 'Color', color);
    plot([cmin{k}(1); cmin{k}(1)], ...
	 [cmin{k}(2); cmax{k}(2)], 'Color', color);
    plot([cmax{k}(1); cmax{k}(1)], ...
	 [cmax{k}(2); cmax{k}(2)], 'Color', color);
    plot([cmin{k}(1); cmax{k}(1)], ...
	 [cmax{k}(2); cmax{k}(2)], 'Color', color);
    plot([cmax{k}(1); cmax{k}(1)], ...
	 [cmin{k}(2); cmax{k}(2)], 'Color', color);
  end
elseif (~strcmp(type, 'none'))
  error(sprintf('Unrecognized plot type: %s', type));
end

% Plot the covariance matrices of all landmarks.
i = map.keySet.iterator;
while (i.hasNext)
  v = i.next;
  p = map.get(v);
  m = p.getMu([]).getArray;
  C = p.getSigma([], []).getArray;
  if (~isempty(R))
    m = R * m(1:2);
    C = R * C(1:2, 1:2) * R';
  end
  if (~isempty(t))
    m = m(1:2) + t;
  end
  % Plot the covariance ellipsoid.
  if (v.equals(f.x))
    color = robot_color;
  else
    id = f.getLandmarkId(v);
    if (~isempty(colors))
      color = colors{id};
    else
      color = 'g';
    end
  end
  plotcov2(m(1:2), C(1:2, 1:2), ...
	   'conf', conf, 'plot-opts', {'Color', color}, ...
	   'num-pts', 20, 'surf-opts', ...
	   {'EdgeAlpha', 0, 'FaceAlpha', 0.1, 'FaceColor', color});
end

