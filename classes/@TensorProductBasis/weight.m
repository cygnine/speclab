function[w] = weight(self, x, d)
% weight -- Evalutes the tensor-product weight function
%
% w = weight(self, x, d)
%
%     Evaluates the weight function for the tensor-product orthogonal
%     polynomial family. The input x should have self.dim columns. The d'th
%     partial derivative is evaluated, where d is an array with self.dim
%     columns. Each row of d specifies a partial derivative:
%
%       d = [3 4 0],   with self.dim = 3
%
%     mean take the 3rd partial wrt the first dimension, the 4th partial wrt
%     the second dimension, and the 0th partial wrt the third dimension.
%
%     The output w is an array with dimensions size(x,1) x size(d,1). Each
%     column of w corresponds to a row of d.

if size(x,2) ~= self.dim
  error('Input points must have the same dimension (columns) as self.dim');
end

if nargin < 3
  d = zeros([1 self.dim]);
elseif size(d,2) ~= self.dim
  error('Specification of partial derivatives must have the same dimension (columns) as self.dim');
end

w = ones([size(x,1) size(d,1)]);
for dd = 1:size(d,1)
  for q = 1:self.dim;
    w(:,dd) = w(:,dd).*self.bases{q}.weight(x(:,q), d(dd,q));
  end
end
