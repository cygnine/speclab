function[x,w] = hermitepoly_gq(N,mu,shift,scale);

% [x,w] = hermitepoly_gq(N,mu,shift,scale);
% Returns the N-point Gaussian quadrature rule for the Hermite polynomials.
% The weight function is \abs{(t-shift)/scale}^{2\mu} e^{-((t-shift)/scale)^2}. 
%
% 20080524: acn

hermite_parameters;

[as,bs] = hermite_recurrence(N+1,mu,shift,scale);

[x,w] = opoly_gq(as,bs,N);
