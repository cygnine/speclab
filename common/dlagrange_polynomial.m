function[l] = dlagrange_polynomial(x, z, n)
% dlagrange_polynomial -- Evaluates the deriative of a Lagrange interpolating polynomial
%
% l = dlagrange_polynomial(x, z, n)
%
%     Evaluates the derivative of the n'th Lagrange interpolating polynomial of
%     the nodes x at the locations z. See lagrange_polynomial for an in-depth
%     explanation of the inputs and outputs, and for examples.

% Again, this implementation loops over n and is therefore slow
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
inds = 1:X;

if Z*X <= maxvarsize  % use lots of storage, is faster
  pix = repmat(z, [1 X]);
  pix = pix - repmat(x.', [Z 1]);

  dpix = zeros([Z 1]);
  for q = 1:X
    tempinds = inds;
    tempinds(q) = [];
    dpix = dpix + prod(pix(:,tempinds),2);
  end

  pix = prod(pix, 2);
else  % forced to use less storage, is *much* slower
  pix = ones(size(z));
  dpix = zeros(size(z));
  for q = 1:X;
    dpix_temp = ones(size(z));
    tempinds = inds;
    tempinds(q) = [];
    for qq = tempinds
      dpix_temp = dpix_temp.*(z-x(qq));
    end
    dpix = dpix + dpix_temp;

    pix = pix.*(z - x(q));
  end
end

l = zeros([Z N]);
ncount = 1;

for nvalue = n
  tempinds = inds([1:nvalue-1 nvalue+1:end]);
  factor = prod(x(nvalue) - x(tempinds));

  l(:,ncount) = 1/factor*(dpix.*(z-x(nvalue)) - pix)./(z-x(nvalue)).^2;

  flags = abs(z-x(nvalue))<tol;

  factor = sum(1./(x(nvalue) - x(tempinds)));
  l(flags,ncount) = factor;

  ncount = ncount+1;
end
