function[P] = primitive_pointwise(x,y, varargin)
% primitive_pointwise -- The pointwise primitive matrix
%
% P = primitive(x, y, {alpha=-1/2, beta=-1/2, normalization='normal', 
%                      shift=0, scale=1})
%
%     Returns a matrix that interpolates point-values specified at x with the
%     Jacobi polynomials of the specified family and returns the primitive
%     of the interpolant evaluated at different points y. The primitive is
%     normalized so that the evaluation at the left-hand endpoint of the
%     Jacobi polynomial interval is 0.

persistent defaults dcoeffs eval_jac sepmat spdiag
if isempty(defaults)
  from speclab.orthopoly.jacobi import defaults
  from speclab.orthopoly.jacobi import jacobi_poly as eval_jac
  from speclab.orthopoly.jacobi import derivative_connection as dcoeffs
  from speclab.orthopoly.jacobi import ...
    integer_separation_connection_matrix as sepmat
  from labtools import spdiag
end

opt = defaults(varargin{:});

N = length(x);
% interpolate to find modal coefficients
V = eval_jac(x, 0:(N-1), opt);
V2 = eval_jac(y, 0:N,opt);

% promote modal coefficients to (alpha+1, beta+1) class
C = sepmat(N, opt.alpha, opt.beta, 1,1);
D = spalloc(N+1,N,nnz(C));

% normalize coefficients to correspond to (alpha,beta) multipliers
etas = dcoeffs(1:N, opt.alpha, opt.beta, opt);
C = spdiag(1./etas)*C;
D(2:end,:) = C;

left_endpoint = opt.shift - opt.scale;

% Find condition to fix left-endpoint
temp = speye(N+1);
temp(1,1) = 0;
temp(1,2:end) = -eval_jac(left_endpoint, 1:N, opt)/eval_jac(left_endpoint,0,opt);

% Yay linear algebra
P = V2*temp*D*inv(V);
