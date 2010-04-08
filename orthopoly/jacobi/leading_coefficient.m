function[k] = leading_coefficient(n,varargin)
% k = leading_coefficient(n, {alpha=-1/2, beta=-1/2, dim=1, normalization='normal', weight_normalization='', shift=0, scale=1})
%
%     Evaluates the leading coefficient of the Jacobi polynomials of specified
%     class with specified normalization.

persistent driver defaults recurrence
if isempty(driver)
  from speclab.orthopoly import leading_poly_coefficient as driver
  from speclab.orthopoly.jacobi import defaults recurrence
end

opt = defaults(varargin{:});

for q = 1:opt.dim
  poly_parameters(q) = struct('alpha', opt.alpha(q), 'beta', opt.beta(q));
end

k = driver(n, recurrence, opt.dim, opt.normalization, opt.weight_normalization, poly_parameters);
