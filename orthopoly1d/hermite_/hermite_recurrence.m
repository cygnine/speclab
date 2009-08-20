function[as,bs] = hermite_recurrence(N,mu,shift,scale)

% [as,bs] = hermite_recurrence(N,mu,shift,scale)
% Calculates the first N recurrence coefficients for the generalized monic 
% hermite polynomials. 
% mu is the polynomial generalization factor: 
% The weight function is \abs{(t-shift)/scale}^{2\mu} e^{-((t-shift)/scale)^2}. 
% 20080524:acn

% Sets default parameter values
hermite_parameters;

as = zeros([N 1]);
bs = as;

bs(1) = gamma(mu+1/2);
oddk = 1:2:(N-1);
bs(2:end) = 1/2*(1:(N-1));
bs(oddk+1) = bs(oddk+1)+mu;
%bs(2:end) = bs(2:end);

% Augments recurrence for affine shifting/scaling
recurrence_scaleshift;
