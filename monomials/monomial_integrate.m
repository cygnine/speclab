function[ints] = monomial_integrate(mc,interval)
% [ints] = monomial_integrate(mc,interval)
%
%     Given modal coefficients (MC) for a monomial expansion, integrates this
%     expansion over the specified interval (INTERVAL). Is vectorized in the
%     columns of mc and interval.

n = size(mc,1);
C = size(mc,2);

ints = zeros([n,C]);

for q = 1:n
  ints(q,:) = 1/q*(interval(2,:).^q - interval(1,:).^q);
end

ints = sum(ints.*mc,1);
