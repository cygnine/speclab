function[a,b] = recurrence(N,varargin)
% recurrence -- recurrence coefficients for Hermite polynomials
%
% [a,b] = recurrence(N, {mu=0})
%     Calculates the first N recurrence coefficients for the generalized hermite
%     polynomials.  

% Sets default parameter values
persistent defaults
if isempty(defaults)
  from speclab.orthopoly1d.hermite import defaults
end

opt = defaults(varargin{:});

a = zeros([N 1]);
b = a;

b(1) = gamma(opt.mu+1/2);
oddk = 1:2:(N-1);
b(2:end) = 1/2*(1:(N-1));
b(oddk+1) = b(oddk+1)+opt.mu;
