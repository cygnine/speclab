function [dd] = divided_difference(x,y)
% [DD] = DIVIDED_DIFFERENCE(X,Y)
%
%     Calculates modal coefficients for interpolation by using the divded
%     difference formulation. The modal coefficients are computed assuming that
%     the basis set is constructed using products of (x-X(k)) for k =
%     1,2,...,n-1 with n being the size of the input. Of course, the basis
%     products are ordered by the input ordering.  The indexing assumed here is
%     along rows. If X and Y have multiple columns (say C of them), this
%     produces a length(X) x C matrix with divided difference coefficients.


[n,C] = size(x);
if and(n==1,C>1)  % I don't think you're calling this for constant interpolants;
  n = C;          % I'm assuming you just transposed some vectors around.
  C = 1;
  x = x';
  y = y';
end

dd = zeros([n,C]);

% First-order differences:
Differences = y;
dd(1,:) = y(1,:);

for q = 2:n
  % straddling x intervals
  RightX = x(q:n,:);
  LeftX = x(1:(n-q+1),:);

  Differences = diff(Differences)./(RightX-LeftX);
  dd(q,:) = Differences(1,:);
end
