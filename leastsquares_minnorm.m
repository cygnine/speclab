function[x] = leastsquares_minnorm(A, b, Cr, Cx)
% leastsquares_minnorm -- Computes the minimum-norm, least-squares solution of a linear system
%
% x = leastsquares_minnorm(A, b, [[ Cr, Cx ]])
%
%     Define the Gram matrices R = Cr^T*Cr and X = Cx^T*Cx. This computes the
%     solution to the optimization problem:
%
%       Among all vectors x minimizing ||A*x - b||_R, find the x minimizing ||x||_X.
%
%     This solution is only unique if the kernels of Cr*A and Cx are
%     non-intersecting. When this solution is not unique, this routine returns
%     the vector x that minimizes the Euclidean norm as well. The solution can
%     be written as a generalized pseudoinverse, and the workhorse for this
%     routine is the gsvd.
%
%     Cr and Cx are the matrix square-root (e.g. Cholesky) factors for the
%     norms R and X. When empty inputs, or no inputs are given, the Euclidean
%     norm is assumed. Supported syntaxes are:
%
%     x = leastsquares_minnorm(A, b)
%     x = leastsquares_minnorm(A, b, Cr, Cx)
%     x = leastsquares_minnorm(A, b, Cr)
%     x = leastsquares_minnorm(A, b, [], Cx)

[m,n] = size(A);

if nargin < 4 || isempty(Cx)
  Cx = speye(n);
  if (nargin < 3) || (isempty(Cr))
    Cr = speye(n);
  end
elseif isempty(Cr)
  Cr = speye(n);
end

% First transform residual so we're doing standard unweighted least-squares:
% A x = b -----> K x = d
d = Cr*b;
K = Cr*A;

% Now get generalized inverse
[u,v,y,c,s] = gsvd(K, Cx);

% Since c only has elements on diagonal max(n-m, 0), just take inverses of that.
cdiag = diag(c, max(n-m,0));
crank = find(cdiag, 1, 'last');
cdiag(1:crank) = 1./cdiag(1:crank);
pinvc = spdiags(cdiag(:), -max(n-m,0), size(c,2), size(c,1)) 

% The generalized inverse:
x = y'\(pinvc*u'*d);
