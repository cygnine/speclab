function[inds] = least_opoly_adaptive_inds(theta, f, N, varargin)
% least_opoly_adaptive_coeffs -- Computes the least orthogonal adaptive interpolant.
%
% c = least_opoly_adaptive_coeffs(theta, f, N, [tol=1e-10, ip=eye(), basis='legendre', initial_inds = 1])
%
%     Blah, blah, explain this.

persistent dim subdim indexing find_order
persistent strict_inputs ip least_opoly_lqlu invl invu
if isempty(dim)
  from labtools import strict_inputs

  from speclab.orthopoly.interp import least_opoly_lqlu 
  from speclab.common.tensor import space_dimension as dim
  from speclab.common.tensor import subspace_dimension as subdim
  from speclab.common.tensor import linear_to_array_indexing as indexing
  from speclab.common.tensor import Npoints_to_poly_order as find_order

  from labtools.linalg import triu_back_substitute as invu
  from labtools.linalg import tril_forward_substitute as invl

  % For default ip, basis:
  ip = @(d,k) speye(subdim(d,k));
end

opt = strict_inputs({'tol', 'basis', 'ip', 'initial_inds'}, {1e-10, [], [], 1}, [], varargin{:});
if isempty(opt.basis)
  bases = cell([size(theta,2) 1]);
  for d = 1:length(bases)
    interval = 1.01*[min(theta(:,d)) max(theta(:,d))];
    bases{d} = LegendrePolynomialBasis('domain', interval);
  end
  TL = TensorProductBasis(bases{:});
  opt.basis = @(x,n) TL.evaluate(x,n);
end
if isempty(opt.ip)
  opt.ip = ip;
end

[Nt,d] = size(theta);

inds = zeros([N 1]);
inds(1:length(opt.initial_inds)) = opt.initial_inds;
all_inds = 1:Nt;
leftover_inds = setdiff(all_inds, inds(1:length(opt.initial_inds)));
for n = (length(opt.initial_inds)+1):N
  possible_inds = leftover_inds;
  vals = zeros(size(possible_inds));
  for m = 1:length(possible_inds)
    newinds = [inds(1:n-1); possible_inds(m)];
    [l2,u2,p2,v2,k2] = least_opoly_lqlu(theta(newinds,:), 'basis', opt.basis);
    a2 = invu(u2, invl(l2, p2*f(newinds)));
    vals(m) = abs(a2(end));
  end
  [garbage, temp] = min(vals);
  inds(n) = possible_inds(temp);
  leftover_inds(temp) = [];
end
%[l,u,p,v,k] = least_opoly_lulq(theta(inds,:), 'basis', opt.basis);

%l = eye(N);
%u = eye(N);
%p = speye(Nt);
%
%interp_inds = zeros([N 1]); % This is a matrix that stores the row indices of
%                            % theta used to interpolate.
%all_inds = (1:Nt).';
%
%% This is just a guess: this vector could be much larger, or much smaller
%v = zeros([1000 1]);
%v_index = 1;
%
%% Current polynomial degree
%k_counter = 0;
%k = zeros([N 1]);  % k(q) gives the degree used to eliminate the q'th point
%
%% The current LU row to factor out:
%lu_row = 1;
%
%%%%%%% Begin k=0 (constant) elimination
%% All initialization is done -- now we have to start with the first point.
%current_dim = 1;
%poly_indices = 1;
%%W = p*opt.basis(theta, poly_indices);
%
%M = opt.ip(d, k_counter);
%sM = chol(M);
%
%% Simple: just find point with largest value
%[garbage, interp_inds(1)] = min(f);
%l(1,1) = sM*opt.basis(theta(interp_inds(1),:), 1);
%%%%%%% End k=0 (constant) elimination
%
%[l,u,p,v,k] = least_opoly_lqlu(
