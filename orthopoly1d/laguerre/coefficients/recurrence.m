function[a,b] = recurrence(N,varargin)
% recurrence -- recurrence coefficients for Laguerre polynomials
%
% [a,b] = recurrence(N, {alpha=0})
%     Calculates the first N recurrence coefficients for the generalized
%     Laguerre polynomials.  

% Sets default parameter values
global handles;
opt = handles.speclab.orthopoly1d.laguerre.defaults(varargin{:});

a = zeros([N 1]);
b = a;

b(1) = gamma(1+opt.alpha);
a = 2*(0:(N-1)).' + opt.alpha + 1;

ks = (1:(N-1)).';
b(2:end) = ks.*(ks+alpha);
