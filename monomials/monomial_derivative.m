function[mc] = monomial_derivative(mc)
% [MC] = MONOMIAL_DERIVATIVE(MC)
%
%     Given the modal coefficients for a polynomial (mc), evaluates the
%     derivative and returns those modal coefficients. Is vectorized in the
%     columns of mc.

n = size(mc,1);

temp = spalloc(n-1,n-1,n-1);
indices = (0:(n-2))*(n-1) + (1:(n-1));
temp(indices) = 1:(n-1);

%mc(1:(n-1),:) = spdiags((1:(n-1)).',0,n-1,n-1)*mc(2:n,:);
mc(1:(n-1),:) = temp*mc(2:n,:);
mc(n,:) = [];

if n==1
  mc = zeros(size(mc));
end
