function[w] = wjacobi_sqrt(r,alpha,beta,shift,scale)

% function[w] = wjacobi_sqrt(r,alpha,beta,shift,scale)
% Implements the square root of the weight function for the Jacobi polynomials at the point x.
% The weight function is (1-(r-shift)/scale)^alpha*(1-(r-shift)/scale)^beta
%
% 20080624: acn

jacobi_parameters;

cr = (r-shift)/scale;
w = (1-cr).^(alpha/2).*(1+cr).^(beta/2);
