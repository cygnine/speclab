% MATLAB File : newton_to_monomial.m
% [mc] = newton_to_monomial(nc,x)
%
% * Creation Date : 2009-06-05
%
% * Last Modified : Fri 12 Jun 2009 03:21:36 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Rewrites the Newton polynomial modal coefficients (nc) as monomial
%   coefficients (mc). Thus, it requires the nodal locations (x) specifying the
%   Newton basis.
%   Is vectorized for multiple columns of nc/x. Assumes size(x) == size(nc),
%   although technically the last row of x is not used. 

function[mc] = newton_to_monomial(nc,x)

n = size(x,1);
C = size(x,2);

mc = zeros(size(nc));
mc = zeros([n,C]);

NewtonBasis = zeros([n,C]);
% First compute monomial contributions from unweighted Newton polys
NewtonBasis(1,:) = 1;
NBTemp = NewtonBasis;

mc(1,:) = nc(1,:);

for k = 2:n
  % Update Newton Basis
  NewtonBasis(1:(k-1),:) = -NBTemp(1:(k-1),:)*spdiags(x(k-1,:).',0,C,C);
  NewtonBasis(2:k,:) = NewtonBasis(2:k,:) + ...
                       NBTemp(1:(k-1),:);
  NBTemp = NewtonBasis;

  % Update monomial coefficients
  mc(1:k,:) = mc(1:k,:) + NewtonBasis(1:k,:)*spdiags(nc(k,:).',0,C,C);
end
