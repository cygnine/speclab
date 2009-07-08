function[w] = wlaguerre(x,alpha,shift,scale)

% function[w] = wlaguerre(x,alpha,shift,scale)
% Implements the weight function for the Laguerre polynomials at the point x.
% The weight function is (t-shift)^alpha*exp(-(t-shift))
%
% 20080623: acn

laguerre_parameters;

cx = (x-shift)/scale;
w = cx.^alpha.*exp(-cx);
