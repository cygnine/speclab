% MATLAB File : monomial_derivative.m
% [mc] = monomial_derivative(mc)
%
% * Creation Date : 2009-06-05
%
% * Last Modified : Fri 12 Jun 2009 03:47:52 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Given the modal coefficients for a polynomial (mc), evaluates the
%   derivative and returns those modal coefficients. Is vectorized in the
%   columns of mc.

function[mc] = monomial_derivative(mc)

n = size(mc,1);

mc(1:(n-1),:) = spdiags((1:(n-1)).',0,n-1,n-1)*mc(2:n,:);
mc(n,:) = [];
