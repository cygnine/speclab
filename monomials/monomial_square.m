% MATLAB File : monomial_square.m
% [mc2] = monomial_square(mc);
%
% * Creation Date : 2009-06-05
%
% * Last Modified : Fri 12 Jun 2009 03:53:23 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Given modal coefficients for a monomial expansion (mc), this
%   function returns the modal coefficients for the function squared, and is
%   vectorized in the columns of mc.
%
%   Note that this functionality should be provided by Matlab's CONV function,
%   but for some strange reason CONV isn't vectorized.

function[mc2] = monomial_square(mc);

n = size(mc,1);
C = size(mc,2);

mc2 = zeros([2*n-1,C]);

for q = 1:n
  mc2((q:(q+n-1)),:) = mc2((q:(q+n-1)),:) + mc*spdiags(mc(q,:).',0,C,C);
end
