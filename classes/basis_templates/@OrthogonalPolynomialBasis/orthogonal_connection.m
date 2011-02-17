function[lmbda] = orthogonal_connection(self,other,N)
% orthogonal_connection -- Connection coefficients between orthogonal families
%
% lambda = orthogonal_connection(self, other, N)
%
%     Given another OrthogonalPolynomialBasis instance 'other', this function
%     uses the recurrence coefficients from both families to compute the
%     connection coefficients for the first N polynomials. Let 'p' be the family
%     'self', and 'q' be the family 'other'. lambda is an 
%     (N x N) lower-triangular matrix with entries that satisfy:
%
%       p_n = \sum_{m=0}^n lambda(n,m) q_m
%
%     where the matrix lambda in the above equation is zero-indexed. The
%     polynomials p and q adopt the normalizations defined by the instances self
%     and other.

persistent spdiag
if isempty(spdiag)
  from labtools import spdiag
end

if not(isa(other, 'OrthogonalPolynomialBasis'))
  error('Input object must be an OrthogonalPolynomialBasis instance');
end

if N < 0
  error('Input N must be a non-negative integer')
elseif N==0
  lmbda = [];
  return
end

N = round(N);

[a,b] = self.recurrence(0:N);
[alpha, beta] = other.recurrence(0:N);

lmbda = zeros(N+1);  % Will chop off extra row + col at end
lmbda(1,1) = sqrt(beta(1)/b(1));
for n = 2:N
  if n==2 % Ignore n=0 term in recurrence
    % For m=1, need to ignore m=0 term
    lmbda(2,1) = 1/sqrt(b(2))*(alpha(1) - a(2))*lmbda(1,1);
    lmbda(2,2) = 1/sqrt(b(2))*((alpha(1) - a(2))*lmbda(1,2) + sqrt(beta(2))*lmbda(1,1));
  else
    % For m=1, need to ignore m=0 term
    lmbda(n,1) = 1/sqrt(b(n))*( -sqrt(b(n-1))*lmbda(n-2,1) + ...
                                (alpha(1) - a(n-1))*lmbda(n-1,1) + ...
                                sqrt(beta(2))*lmbda(n-1,2) );

    % Can vectorize the rest
    lmbda(n,2:n) = 1/sqrt(b(n))*(-sqrt(b(n-1))*lmbda(n-2,2:n) + ...
                                  sqrt(beta(2:n)).*lmbda(n-1,1:n-1) + ...
                                  (alpha(2:n) - a(n-1)).*lmbda(n-1,2:n) + ...
                                  sqrt(beta(3:n+1)).*lmbda(n-1,3:n+1));
  end
end

lmbda = lmbda(1:N,1:N);

% Rescale according to normalization
k = self.scale_functions(ones([1 N]), 0:(N-1));
kappa = other.scale_functions(ones([1 N]), 0:(N-1));
lmbda = spdiag(k)*lmbda*spdiag(1./kappa);

