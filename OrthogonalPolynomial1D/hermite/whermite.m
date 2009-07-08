function[w] = whermite(x,mu,shift,scale)

% function[w] = whermite(x,mu,shift,scale)
% Implements the weight function for the Hermite polynomials at the point x.
% The weight function is \abs{(t-shift)/scale}^{2\mu} e^{-((t-shift)/scale)^2}. 
%
% 20080624: acn

hermite_parameters;

cx = 1/scale*(x-shift);
w = abs(cx).^(2*mu).*exp(-cx.^2);
