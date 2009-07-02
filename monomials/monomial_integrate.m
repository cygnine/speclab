% MATLAB File : monomial_integrate.m
% [ints] = monomial_integrate(mc,interval)
%
% * Creation Date : 2009-06-05
%
% * Last Modified : Fri 12 Jun 2009 03:52:30 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Given modal coefficients (mc) for a monomial expansion, integrates this
%   expansion over the specified interval (interval). Is vectorized in the
%   columns of mc and interval.

function[ints] = monomial_integrate(mc,interval)

n = size(mc,1);
C = size(mc,2);

ints = zeros([n,C]);

for q = 1:n
  ints(q,:) = 1/q*(interval(2,:).^q - interval(1,:).^q);
end

ints = sum(ints.*mc,1);
