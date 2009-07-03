function[mc2] = monomial_square(mc);
% [MC2] = MONOMIAL_SQUARE(MC);
%
%     Given modal coefficients for a monomial expansion (mc), this
%     function returns the modal coefficients for the function squared, and is
%     vectorized in the columns of mc.
%
%     Note that this functionality should be provided by Matlab's CONV function,
%     but for some strange reason CONV isn't vectorized.
%
%     For vectors:
%     n = length(MC), and is the (poly order + 1)

n = size(mc,1);
C = size(mc,2);

mc2 = zeros([2*n-1,C]);

for q = 1:n
  mc2((q:(q+n-1)),:) = mc2((q:(q+n-1)),:) + mc*spdiags(mc(q,:).',0,C,C);
end
