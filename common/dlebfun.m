function[y] = dlebfun(x,z)
% dlebfun -- Evaluates the derivative of the Lebesgue function
%
% y = dlebfun(x,z)
%
%     For a vector x of nodal locations, this evaluates the derivative of the
%     associated Lebesgue function at the locations z. See lebfun for a more
%     detailed description and examples.

% This implementation loops over all `intervals', over each of which the
% function is a polynomial.

persistent dlagrange_polynomial
if isempty(dlagrange_polynomial)
  from speclab.common import dlagrange_polynomial
end

x = sort(x(:)); X = length(x);
zsize = size(z);
z = z(:);

L = dlagrange_polynomial(x, z, 1:X);

[garbage, bin] = histc(z, [-Inf; x; Inf]);
bin(bin==X+2) = X+1;

inds = 1:X;
signs_rightplus = (-1).^(inds-X);
signs_leftplus = (-1).^(inds-1);

y = zeros([length(z) 1]);
for q = 0:X;  % Loop over each interval ( x(q), x(q+1) )
  flags = bin==(q+1);

  if q==0
    signs = signs_leftplus;
  elseif q==X
    signs = signs_rightplus;
  else
    signs = zeros([1 X]);
    len = q;
    signs(1:q) = signs_rightplus(end-q+1:end);
    len = X-q;
    signs((q+1):X) = signs_leftplus(1:len);
  end

  y(flags) = L(flags,:)*signs(:);
end

y = reshape(y, zsize);
