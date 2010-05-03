function[l] = lagrange_polynomial(x, z, n)
% lagrange_polynomial -- Evaluates a Lagrange interpolating polynomial of a nodal set
%
% l = lagrange_polynomial(x, z, n)
%
%     x is a collection of one-dimensional nodes of length N. n is an indexing
%     vector containing values between 1 and N. The matrix returned l has size
%     length(z) x length(n) and each column evaluates the n'th Lagrange
%     polynomial interpolant from nodal set x at the locations z.
%
%     If n is a scalar, then the output l has the same size as z.
%
%
%     Evaluate the fourth Lagrange interpolant of 7 equidistant nodes:
%      x = linspace(-1, 1, 7);
%      z = linspace(-1, 1, 1e3);
%      l = lagrange_polynomial(x, z, 4);
%      plot(z, l, 'k', x(4), 1, 'r.', x([1:3 5:end]), 0, 'r.');
%
%     Evaluates three of the Lagrange interpolants of the same nodes:
%      l = lagrange_polynomial(x, z, [1 3 6]);
%      plot(z, l, '-', x([1 3 6]), 1, 'k.', x([2 4 5 7:end]), 0, 'k.')

% The implementation loops over n
n = n(:).';  N = length(n);
x = sort(x(:)); X = length(x);
zsize = size(z);
z = z(:); Z = length(z);

% Error checking:
if any(n>X) | any(n<1)
  error('The interpolant indices must be values 1, 2, ..., length(x)');
end

maxvarsize = 1e8;
tol = 1e-13;

if Z*X <= maxvarsize  % use lots of storage, is faster
  pix = repmat(z, [1 X]);
  pix = pix - repmat(x.', [Z 1]);
  pix = prod(pix, 2);
else  % forced to use less storage, is slower
  pix = ones(size(z));
  for q = 1:X;
    pix = pix.*(z - x(q));
  end
end

l = repmat(pix, [1 N]);
inds = 1:X;

% This is slow if N is large
ncount = 1;
for nvalue = n;
  tempinds = inds([1:nvalue-1 nvalue+1:end]);
  factor = prod(x(nvalue) - x(tempinds));
  l(:,ncount) = l(:,ncount)./(factor*(z-x(nvalue)));

  flags = abs(z-x(nvalue))<tol;
  l(flags,ncount) = 1;
  ncount = ncount+1;
end

if length(n)==1
  l = reshape(l, zsize);
end
