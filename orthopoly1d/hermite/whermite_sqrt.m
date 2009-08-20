function[w] = whermite_sqrt(x,mu,shift,scale)

% function[w] = whermite_sqrt(x,mu,shift,scale)
% Implements the square root of the weight function for the Hermite polynomials at the point x.
% This square root is a non-principal branch cut. 
% The weight function is \abs{(t-shift)/scale}^{2\mu} e^{-((t-shift)/scale)^2}. 
%
% 20080624: acn

hermite_parameters;

cx = 1/scale*(x-shift);
w = cx.^mu.*exp(-1/2*cx.^2);
