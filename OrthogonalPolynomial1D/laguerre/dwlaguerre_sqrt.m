function[w] = dwlaguerre_sqrt(x,alpha,shift,scale)

% function[w] = dwlaguerre_sqrt(x,alpha,shift,scale)
% Implements the derivative of the square root of the weight function for the Laguerre polynomials at the point x.
% The weight function is ((t-shift)/scale)^alpha*exp(-(t-shift)/scale)
%
% 20080623: acn

laguerre_parameters;

cx = (x-shift)/scale;

w = wlaguerre_sqrt(x,alpha-2,shift,scale);
w = w.*(alpha-cx)/(2*scale);
