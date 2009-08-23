function[x,w] = gauss_lobatto_quadrature(N,varargin)
% [x,w] = gauss_lobatto_quadrature(N, {alpha=-1/2,beta=-1/2,shift=0,scale=1,r1=-1,r2=1})
%
%     Returns the N-point Gauss-Lobatto quadrature rule for the Jacobi polynomials.
%     The weight function is given by speclab.orthopoly1d.jacobi.weights.weight.
%     The Lobatto points are located at x=r1,r2.
%
%     TODO: Chebyshev case exception

global handles;
opoly = handles.speclab.orthopoly1d;
jac = opoly.jacobi;
pss = handles.speclab.common.physical_scaleshift_1d;
sss = handles.speclab.common.standard_scaleshift_1d;

opt = jac.defaults(varargin{:});
[alpha,beta,scale,shift,r1,r2] = ...
  deal(opt.alpha,opt.beta,opt.scale,opt.shift,opt.r1,opt.r2);

% Move all inputs to the standard interval [-1,1]
r1 = sss(r1,opt);
r2 = sss(r2,opt);

% Compute recurrence constants
[a,b] = jac.coefficients.recurrence(N,opt);

% Solve eigenvalue problem
[x,w] = opoly.gauss_lobatto_quadrature(a,b,N,r1,r2);

% Convert back to physical inteval
x = pss(x,opt);
w = w*scale;
