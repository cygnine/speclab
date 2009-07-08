function[x,w] = gauss_lobatto_quadrature(N,varargin)
% [X,W] = GAUSS_LOBATTO_QUADRATURE(N,
%             {ALPHA=-1/2,BETA=-1/2,SHIFT=0,SCALE=1,R1=-1,R2=1})
%
%     Returns the N-point Gauss-Lobatto quadrature rule for the Jacobi polynomials.
%     The weight function is
%     (1-1/SCALE*(x-SHIFT))^ALPHA*(1+1/SCALE*(x-SHIFT))^BETA. The Lobatto points
%     are located at x=R1,R2.
%
%     TODO: Chebyshev case exception

global handles;
opoly = handles.speclab.OrthogonalPolynomial1D;
jac = opoly.jacobi;
pss = handles.speclab.common.physical_scaleshift_1d;
sss = handles.speclab.common.standard_scaleshift_1d;

jac.defaults;
[alpha,beta,scale,shift,r1,r2] = ...
  deal(opt.alpha,opt.beta,opt.scale,opt.shift,opt.r1,opt.r2);

% Move all inputs to the standard interval [-1,1]
r1 = sss(r1,scale,shift);
r2 = sss(r2,scale,shift);

% Compute recurrence constants
[a,b] = jac.recurrence(N,alpha,beta,shift,scale);

% Solve eigenvalue problem
[x,w] = opoly.gauss_lobatto_quadrature(a,b,N,r1,r2);

% Convert back to physical inteval
x = pss(x,scale,shift);
w = w*scale;
