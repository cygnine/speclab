function[f] = clenshaw_algorithm(x, coeffs, varargin)
% clenshaw_algorithm -- Evaluates a Jacobi polynomial sum
%
% f = clenshaw_algorithm(x, coeffs, {alpha=-1/2, beta=-1/2, normalization='normal', scale=1, shift=0})
%
%     Uses the Clenshaw algorithm to evaluate the sum 
%
%       f(x) = \sum_n coeffs(n)*p_{n-1}(x),
%
%     where p_n are Jacobi polynomials of class (alpha,beta) with specified
%     normalization ('normal' or 'monic').

persistent clenshaw recurrence defaults
if isempty(recurrence)
  from speclab.orthopoly.jacobi import recurrence defaults
  from speclab.orthopoly import clenshaw_evaluation as clenshaw
end

opt = defaults(varargin{:});
[alpha,beta] = recurrence(length(coeffs)+1, opt);

f = clenshaw(x, coeffs, alpha, beta, opt);
