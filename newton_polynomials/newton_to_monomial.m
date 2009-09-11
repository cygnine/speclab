function[mc] = newton_to_monomial(nc,x)
% [MC] = NEWTON_TO_MONOMIAL(NC,X)
%
%     Rewrites the Newton polynomial modal coefficients (NC) as monomial
%     coefficients (MC). Thus, it requires the nodal locations (X) specifying
%     the Newton basis.  Unfortunately, the standardization I chose regarding
%     order of the monomial coefficients is opposite to that of Matlab's
%     POLYVAL, so you'll have to flipud things if you want to use that function
%     to transform polynomials. 
%
%     Is vectorized for multiple columns of NC/X. Assumes size(X) == size(NC),
%     although technically the last row of X is not used. 

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
