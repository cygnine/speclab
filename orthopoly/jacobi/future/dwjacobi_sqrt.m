function[w] = dwjacobi_sqrt(x,alpha,beta,shift,scale)

% function[w] = dwjacobi_sqrt(x,alpha,beta,shift,scale)
% Implements the derivative of the square root of the weight function for the Jacobi polynomials at the point x.
% The weight function is (1-(r-shift)/scale)^alpha*(1-(r-shift)/scale)^beta
%
% 20080624: acn

jacobi_parameters;

cx = (x-shift)/scale;
w = (-alpha/2*(1+r) + beta/2*(1-r)).*wjacobi(x,(alpha/2)-1,(beta/2)-1,shift,scale)/scale;
