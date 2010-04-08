function[C] = inv_monomial_connection(N, varargin)
% inv_monomial_connection -- Computes the monomial-Jacobi connection matrix
%
% C = inv_monomial_connection(N, recurrence, {dim=1, alpha=-1/2, beta=-1/2, normalization='normal', weight_normalization=''})
%
%     Computes the N x N lower triangular connection matrix between monomials
%     and Jacobi polynomials, i.e.
%
%     p_{n-1} = \sum_{m=1}^{n} C(n,m) x^{m-1}
%
%     where p_m is the m'th degree Jacobi polynomial. 

persistent defaults recurrence driver
if isempty(defaults)
  from speclab.orthopoly.jacobi import recurrence defaults
  from speclab.orthopoly import inv_monomial_connection_driver as driver
end

opt = defaults(varargin{:});
recurrences = {};

for d = 1:opt.dim;
  recurrences{d} = @(n) recurrence(n, 'alpha', opt.alpha(d), 'beta', opt.beta(d));
end

dopt.dim = opt.dim;
dopt.normalization = opt.normalization;
dopt.weight_normalization = opt.weight_normalization;

C = driver(N, recurrences, dopt);
