function[x,w] = clenshaw_curtis(N)
% clenshaw_curtis -- Returns a Clenshaw-Curtis quadrature rule
%
% [x,w] = clenshaw_curtis(N)
%
%     Returns the N-point clenshaw curtis quadrature rule for integration on the
%     interval [-1,1].
%
%        \int_{-1}^1 f \dx{x} \approx  sum(w.*f(x))

% The nodes are given by the size-N Chebyshev Gauss-Lobatto points. 
x = cos(((N-1):(-1):0).'*pi/(N-1));

% The nodes arccos(x) are nodes for a type-I DCT.

% First generate the quadrature matrix that takes the evaluations at x to the
% cosine coefficients.

[ks, ns] = ndgrid(0:2:(N-1), 0:(N-1));

ak = zeros(N);

ak = cos(ks.*ns*pi/(N-1));
ak(:,1) = 1/2*ak(:,1);
ak(:,end) = 1/2*ak(:,end);
ak = ak*(2/(N-1));

% And these are the weights for each coefficient:
ks = 0:2:(N-1);
ak_weights = ones(size(ks));
if mod(N,2)==0
  ak_weights(2:end) = 2./(1 - (ks(2:end)).^2);
else
  ak_weights(2:end-1) = 2./(1 - (ks(2:end-1)).^2);
  ak_weights(end) = 1/(1 - (N-1)^2);
end

% Put things together:
w = (ak_weights*ak).';
