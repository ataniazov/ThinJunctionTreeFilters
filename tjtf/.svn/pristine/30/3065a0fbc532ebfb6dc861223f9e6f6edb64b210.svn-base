% REGULARIZE - Regularizes a matrix so that it is positive definite.
%
%
%
%
function s = regularize(s)

if (isempty(s))
  return;
end
thresh = sqrt(eps);
span = 1e5;
[V, D] = eig(full(s));
% Eliminate imaginary and negative eigenvalues.
d = real(diag(D));
d(find(d < 0)) = 0;
% Compute the condition number of the matrix.
if ~min(d)
  cond = Inf;
else
  cond = max(d) / min(d);
end
% If the condition number is greater than a threshold, then
% regularize the matrix by filling out its eigenspectrum.
if (cond > span)
  min_eig = max(d) / span;
  d = d + min_eig;
end
s = V * diag(d) * V';
