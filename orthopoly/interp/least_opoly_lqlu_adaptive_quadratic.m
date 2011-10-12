function[l,u,p,v,k,inds,cf] = least_opoly_lqlu_adaptive_quadratic(theta, f, N, varargin)
% least_opoly_lqlu_adpative_quadratic -- Computes de Boor's LU factorization with LQ factorizations
%
% [l,u,v,p,k,inds,new_f] = least_opoly_lqlu_adaptive_quadratic(theta, f, N, [tol=1e-10, ip=eye(),  basis='legendre', initial_inds = 1])
%
%     Using the sequential method, this routine chooses N points using the
%     (theta,f) data that minimizes the change in the coefficients. 
%
%     This function first computes a least-squares quadratic, and then uses that
%     to start the sequential addition process.
%
%     We always start out with an initial collection of nodes. The matrix
%     initial_inds determines this.

persistent dim strict_inputs subdim ip eyecols
persistent least_opoly_lqlu_sequential_adaptive
if isempty(dim)
  from speclab.common.tensor import space_dimension as dim
  from speclab.common.tensor import subspace_dimension as subdim
  from labtools import strict_inputs eyecols
  from speclab.orthopoly.interp import least_opoly_lqlu_sequential_adaptive

  ip = @(d,k) speye(subdim(d,k));
end

opt = strict_inputs({'tol', 'basis', 'ip'}, {1e-10, [], []}, [], varargin{:});
if isempty(opt.basis)
  bases = cell([size(theta,2) 1]);
  for d = 1:length(bases)
    interval = [min(theta(:,d)) max(theta(:,d))];
    bases{d} = LegendrePolynomialBasis('domain', interval);
  end
  TL = TensorProductBasis(bases{:});
  opt.basis = @(x,n) TL.evaluate(x,n);
end
if isempty(opt.ip)
  opt.ip = ip;
end

d = size(theta, 2);

% First compute least-squares quadratic
V = opt.basis(theta, 1:dim(d,1));
c_lsq = V\f; % least-squares solution

cf = [c_lsq; f];
Nq = length(c_lsq);
cx = [zeros([Nq d]); theta];

newbasis = @(y,n) [eyecols(Nq,n); opt.basis(y(Nq+1:end,:), n)];

copt = opt;
copt.basis = newbasis;
copt.initial_inds = 1:Nq;

[l,u,p,v,k,inds] = least_opoly_lqlu_sequential_adaptive(cx, cf, N+Nq, copt);
