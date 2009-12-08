function[inds] = poly_order_to_linear(n, varargin)
% poly_order_to_linear -- transforms a scalar maximum order to linear indexing
%
% inds = poly_order_to_linear(n, {dim=1})
%
%     Given a maximal polynomial order expansion on a tensorized basis, returns
%     all the linear indices necessary to have an expansion up to that order.
%     The dimension of the space is given by dim.

persistent strict_inputs
if isempty(strict_inputs)
  from labtools import strict_inputs
end

opt = strict_inputs({'dim'}, {1}, [], varargin{:});

maxn = nchoosek(n+opt.dim, opt.dim);

inds = (0:(maxn-1)).';
