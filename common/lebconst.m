function[L] = lebconst(x, varargin)
% lebconst -- Computes the polynomial Lebesgue constant of one-dimensional nodes
%
% L = lebconst(x, {[-1,1]})
%
%     Given N nodes x in one dimension this function computes the Lebesuge
%     constant for polynomial interpolation. The Lebesgue constant is computed,
%     by default, for the interval [-1,1]. A different interval may be specified
%     by supplying an optional argument.
%
%     
%     Compute the Lebesgue constant of 7 equidistant nodes;
%      x = linspace(-1, 1, 7);
%      L = lebconst(x);
%
%     Compute the Lebesgue constant for the same nodes, but a different interval:
%      L = lebconst(x, [-2,2]);

persistent bisection dlebfun lebfun
if isempty(bisection)
  from labtools.rootfind import bisection
  from speclab.common import dlebfun lebfun
end

x = sort(x(:)); X = length(x);
f = @(z) dlebfun(x,z);

% Define the interval of approximation
if nargin>1
  interval = varargin{1};
else
  interval = [-1, 1];
end

% Allocation locations where Lebesgue function has local maxima
L_argmaxes = zeros([X+1 1]);

interval_lengths = diff(x);

% Starting locations for root-finding
x0s = x(1:(X-1)) + 1e-8*interval_lengths;
x1s = x(2:X) - 1e-8*interval_lengths;

% Locations of maximum of function g are where g' vanishes
L_argmaxes(1:(X-1)) = bisection(x0s, x1s, f);

% Add endpoints of interval to find maxima
L_argmaxes(X:end) = interval(:);

% Evaluate Lebesgue function at these location and return maximum
L = max(lebfun(x, L_argmaxes));
