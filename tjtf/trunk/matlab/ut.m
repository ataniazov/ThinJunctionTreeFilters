%  UT - The Unscented Transform method of approximating noisy nonlinear
%       functions of Gaussians as affine with Gaussian noise.  
%
%       Let f(x, z) be a vector-valued function where z is a
%       white-noise vector variable.  Then this function computes
%       an approximation of the form
%   
%                    y = g(x) = a + B * x + w,   w ~ N(0, G)
%
%       where a, B, and G are the parameters.
%
%       The approximation is computed as follows. A set of points that
%       are representative of the input's mean and covariance (called
%       sigma points) are chosen, passed through the nonlinearity, and
%       then used to estimate the conditional mean and covariance of the
%       distribution at the other side.  (See the reference below for
%       details.)  Then this conditional distribution is translated
%       into the representation of g(x) above.
%
%       This function counts flops.  (See FLOPS.)
%
% Usage:
% 
%   [a, B, G, ym, yC, CC, sp] = ut(m, C, Q, f[, OPTIONS])
%
% Inputs:
%
%   m    - an n x 1 vector giving the mean of the input vector
%   C    - an n x n symmetric positive semi-definite matrix giving
%          the covariance of the input vector (this 
%          can be [] to indicate no noise).
%   Q    - a k x k symmetric positive semi-definite matrix giving
%          the covariance of the input white noise vector z (this 
%          can be [] to indicate no noise).
%   f    - a function handle to a non-linear function of the form
%                             y = f(x, z, ...)
%          given input x and white noise vector z.  f is assumed to be
%          vectorized so that x and z can be matrices whose columns
%          are individual values (in which case its output must be
%          a matrix whose columns are individual values).
%
% Options:
%
%   'args'    - a cell array of auxiliary arguments that are passed
%               to f.
%   'alpha'
%   'beta' 
%   'kappa'   - These are three parameters used to select sigma points;
%               their default values are optimal for truly Gaussian
%               distributed state.  (alpha is a scaling parameter,
%               beta is the extra weight on the first sigma point, and
%               kappa is a second scaling parameter.  See the provided
%               reference for details on these parameters.)
%
% Outputs:
%   a   - an n x 1 vector giving the constant term of g(x)
%   B   - an n x m symmetric positive semi-definite matrix giving
%         the linear term of g(x)
%   G   - a n x k symmetric positive semi-definite matrix giving
%         the linear noise model of g(x)
%   ym  - the mean of y
%   yC  - the covariance of y
%   CC  - the cross-covariance between y and x
%   sp  - the joint (x, z) sigma points
%
% Reference: Simon J. Julier and Jeffrey K. Uhlmann, "A New Extension
%            of the Kalman Filter to Nonlinear Systems."  In
%            Proceedings of AeroSense: The 11th International
%            Symposium on Aerospace/Defense, Simulation and Controls,
%            Orlando, Florida, 1997.
%
% See also:  UKF
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

function [a, B, G, ym, yC, CC, usp] = ut(m, C, Q, f, varargin)

import javaslam.util.Flops;

n = length(m);   % size of input vector
k = size(Q, 2);  % size of noise vector

no_noise = 0;
if (isempty(C))
  no_noise = 1;
  C = zeros(length(m), length(m));
end

C = regularize(C);
Q = regularize(Q);

[alpha, ...
 beta, ...
 kappa, ...
 args] = process_options(varargin, 'alpha', 1, ...
			           'beta', 0, ...
			           'kappa', 0, ... %3 - n, ...
			           'args', {});

% Augment the state space to include the noise vector; call this new
% state space U.  By setting the covariances between these vectors to
% be zero, we ensure they are independent (as is required by the
% model).
um = [m(:); zeros(k, 1)];
uC = [C,          zeros(n, k); ...
      zeros(k, n),     Q];

% Compute the sigma points (and their weights) for this augmented
% distribution.
[usp, uspmw, uspcw] = sigma_pts(um, uC, alpha, beta, kappa);

% Send the sigma points through the non-linearity
ysp = feval(f, usp(1:n, :), usp((n+1):(n+k), :), args{:});

% Compute the mean and covariance of the images of the sigma points.
[ym, yC] = wsp2mc(ysp, uspmw, uspcw);

% Using the weighted sigma points, compute an estimate of the cross
% covariance between the input and output.
CC = wsp2cc(usp(1:n, :), ysp, uspmw, uspcw);

% Reparameterize.  We know CC = C * B', so B = (C \ CC)'; also,
% yC = B * C * B' + G, so G = yC - B * C * B' = yC - B * CC.
if no_noise
  B = zeros(size(CC, 2), size(C, 1));
else
  B = (C \ CC)';
end
a = ym - B * m;
G = yC - B * CC;

% Count flops.
if (~no_noise)
  Flops.count(Flops.solve(size(C, 1), size(C, 2), size(CC, 2), 1));
end;
Flops.count(length(ym) + ...
	    Flops.mult(size(B, 1), size(B, 2), 1) + prod(size(yC)) + ...
	    Flops.mult(size(B, 1), size(B, 2), size(CC, 2)));

% Regularize G if necessary.
G = regularize(G);

return

% Computes a set of weighted sigma points to represent a multivariate
% Gaussian distribution.  m is a column vector representing the mean,
% C is a covariance matrix, and alpha, beta and kappa are scaling
% parameters.  sp is a matrix whose columns are sigma points, and wm
% and wc are vectors of weights that can be used to compute the mean
% and covariance from the sigma points.
function [sp, wm, wc] = sigma_pts(m, C, alpha, beta, kappa)
  import javaslam.util.Flops;
  n = length(m);    % Dimension of variable.
  npts = 2 * n + 1; % Number of sigma points.
  lambda = alpha^2 * (n + kappa) - n;  % a scaling parameter
  % Calculate matrix square root of weighted covariance matrix.
  C = sparse(C);
  if (issparse(C))
    Cr = cholinc((n + lambda) * C, 'inf')';
    Cr(find(isinf(Cr))) = 0;
  else
    Cr = chol((n + lambda) * C)';  
  end
  % Create the array of the sigma points.
  tmp = repmat(m, [1 n]);
  sp = [m, (tmp - Cr), (tmp + Cr)];
  % Create the array of the weights.
  wm = [lambda, 0.5 * ones(1, npts - 1)] / (n + lambda);
  wc = [(wm(1) + 1 - alpha^2 + beta), wm(2:end)];
  % Count flops.
  Flops.count(11 + 3 * prod(size(C)) + Flops.chol(size(C, 1)));
  return;

% Computes the mean and covariance of a multivariate Gaussian
% distribution identified by a set of weighted sigma points.  The
% rows of sp are the sigma points, and w is a column vector.
function [M, C] = wsp2mc(sp, wm, wc)
  import javaslam.util.Flops;
  [m, n] = size(sp);
  % Calculate the weighted average of the sigma points.
  M = sum(sp .* repmat(wm, [m 1]), 2);
  Flops.count(2 * prod(size(sp)) + m); % Count flops.
  % Center the sigma points.
  X = sp - repmat(M, [1 n]);
  Flops.count(prod(size(sp))); % Count flops.
  % Compute the covariance
  C = sparse((repmat(wc, [m 1]) .* X) * X');
  Flops.count(prod(size(X)) + Flops.mult(m, n, m)); % Count flops.
  return

% Computes the cross-covariance matrix of two sets of sigma points
% with one set of weights.  The rows of sp1 and sp2 are the sigma
% points, and w is a column vector of weights.
function CC = wsp2cc(sp1, sp2, wm, wc)
  import javaslam.util.Flops;
  [u, n] = size(sp1); 
  [v, n] = size(sp2); 
  % Calculate the weighted average of the sigma points.
  M1 = sum(sp1 .* repmat(wm, [u 1]), 2);
  M2 = sum(sp2 .* repmat(wm, [v 1]), 2);
  Flops.count(2 * (prod(size(sp1)) + prod(size(sp2)))); % Count flops.
  % Center the sigma points.
  X1 = sp1 - repmat(M1, [1 n]);
  X2 = sp2 - repmat(M2, [1 n]);
  Flops.count(prod(size(sp1)) + prod(size(sp2))); % Count flops.
  % Compute the cross covariance
  CC = sparse((repmat(wc, [u 1]) .* X1) * X2');
  Flops.count(prod(size(X1)) + Flops.mult(u, n, v)); % Count flops.
  return
