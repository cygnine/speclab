function[x,w] = gauss_lobatto_quadrature(N,varargin)
% [x,w] = gauss_lobatto_quadrature(N, {alpha=-1/2,beta=-1/2,shift=0,scale=1,r1=-1,r2=1})
%
%     Returns the N-point Gauss-Lobatto quadrature rule for the Jacobi polynomials.
%     The weight function is given by speclab.orthopoly.jacobi.weights.weight.
%     The Lobatto points are located at x=r1,r2.
%
%     TODO: Chebyshev case exception

persistent defaults sss pss glq recurrence
if isempty(defaults)
  from speclab.orthopoly.jacobi import defaults
  from speclab.orthopoly.jacobi import recurrence
  from speclab.common import physical_scaleshift_1d as pss
  from speclab.common import standard_scaleshift_1d as sss
  from speclab.orthopoly import gauss_lobatto_quadrature as glq
end

opt = defaults(varargin{:});
[alpha,beta,scale,shift,r1,r2] = ...
  deal(opt.alpha,opt.beta,opt.scale,opt.shift,opt.r1,opt.r2);

% Move all inputs to the standard interval [-1,1]
r1 = sss(r1,opt);
r2 = sss(r2,opt);

% Compute recurrence constants
[a,b] = recurrence(N,opt);

% Solve eigenvalue problem
[x,w] = glq(a,b,N,r1,r2);

% Convert back to physical inteval
x = pss(x,opt);
w = w*scale;
