function[x,w] = gauss_lobatto_quadrature(N,varargin)
% [x,w] = gauss_lobatto_quadrature(N, {alpha=-1/2,beta=-1/2,shift=0,scale=1,r1=-1,r2=1})
%
%     Returns the N-point Gauss-Lobatto quadrature rule for the Jacobi polynomials.
%     The weight function is given by speclab.orthopoly.jacobi.weights.weight.
%     The Lobatto points are located at x=r1,r2.
%
%     TODO: Chebyshev case exception

persistent defaults sss pss glq recurrence gq tensorize_vectors tensorize
if isempty(defaults)
  from speclab.orthopoly.jacobi import defaults
  from speclab.orthopoly.jacobi import recurrence
  from speclab.common import physical_scaleshift as pss
  from speclab.common import standard_scaleshift_1d as sss
  from speclab.orthopoly import gauss_lobatto_quadrature as glq
  from speclab.orthopoly import gauss_quadrature as gq
  from speclab.common.tensor import tensorize_vectors tensorize
end

opt = defaults(varargin{:});
[alpha,beta,scale,shift,r1,r2] = ...
  deal(opt.alpha,opt.beta,opt.scale,opt.shift,opt.r1,opt.r2);

% Move all inputs to the standard interval [-1,1]
r1 = sss(r1,opt);
r2 = sss(r2,opt);

% Compute recurrence constants
xs = cell([opt.dim 1]);
ws = cell([opt.dim 1]);
for q = 1:opt.dim
  [a,b] = recurrence(N+1, 'alpha', opt.alpha(q), 'beta', opt.beta(q));

  % Solve eigenvalue problem
  if N>1
    [xs{q},ws{q}] = glq(a,b,N,r1(q,1),r2(q,1));
  else
    [x{q},ws{q}] = gq(a,b,N);
  end
end
% 1D version:
%[a,b] = recurrence(N,opt);

if opt.dim>1
  x = tensorize_vectors(xs{:});
  w = prod(tensorize_vectors(ws{:}), 2);
else
  x = xs{1};
  w = ws{1};
end


% Convert back to physical inteval
x = pss(x,opt);
w = w*prod(scale);
