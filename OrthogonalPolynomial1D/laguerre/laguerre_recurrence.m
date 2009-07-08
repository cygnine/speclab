function[as,bs] = laguerre_recurrence(N,alpha,shift,scale)

% [as,bs] = laguerre_recurrence(N,alpha,shift,scale)
% Calculates the first N recurrence coefficients for the generalized monic 
% Laguerre polynomials. alpha is the generalized parameter (alpha>-1)
% The weight function is (1/scale*(t-shift))^{alpha}*exp(-1/scale*(t-shift)), the interval is [shift, \infty)
% 20080623:acn

% Sets default parameter values
laguerre_parameters;

as = zeros([N 1]);
bs = as;

%as = C*as;

bs(1) = gamma(1+alpha);

as = 2*(0:(N-1)).'+ alpha + 1;

ks = (1:(N-1)).';
bs(2:end) = ks.*(ks+alpha);

% Take care of shift, scale
recurrence_scaleshift;
