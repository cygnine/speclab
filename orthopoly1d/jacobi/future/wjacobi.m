function[w] = wjacobi(r,alpha,beta,shift,scale)

% function[w] = wjacobi(r,alpha,beta,shift,scale)
% Implements the weight function for the Jacobi polynomials at the point x.
% The weight function is (1-(r-shift)/scale)^alpha*(1-(r-shift)/scale)^beta
%
% 20080624: acn

jacobi_parameters;

cx = (r-shift)/scale;
w = (1-r).^alpha.*(1+r).^beta;
