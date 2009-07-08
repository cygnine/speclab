function[w] = wlaguerre_sqrt(x,alpha,shift,scale)

% function[w] = wlaguerre_sqrt(x,alpha,shift,scale)
% Implements the square root of the weight function for the Laguerre polynomials at the point x.
% The weight function is ((t-shift)/scale)^alpha*exp(-(t-shift)/scale)
%
% 20080623: acn

laguerre_parameters;

cx = (x-shift)/scale;
w = cx.^(alpha/2).*exp(-cx/2);
