function[x,w] = hermitefun_pgq(N,mu,shift,scale);

% [x,w] = hermitefun_pgq(N,mu,shift,scale);
% Returns the N-point pi-Gaussian quadrature rule for the Hermite functions,
% orthogonal over the real line with unit weight.
% The weight function is \abs{(t-shift)/scale}^{2\mu} e^{-((t-shift)/scale)^2}. 
%
% 20080524: acn

hermite_parameters;

[as,bs] = hermite_recurrence(N+1,mu,shift,scale);

[x,w] = opoly_gq(as,bs,N);

w = w./whermite(x,mu,shift,scale);
