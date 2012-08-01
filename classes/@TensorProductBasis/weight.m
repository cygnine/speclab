function[w] = weight(self, x)
% weight -- Evalutes the tensor-product weight function
%
% w = weight(self, x)
%
%     Evaluates the weight function for the tensor-product orthogonal
%     polynomial family. The input x should have self.dim columns.

if size(x,2) ~= self.dim
  error('Input points must have the same dimension (columns) as self.dim');
end

w = ones([size(x,1) 1]);
for q = 1:self.dim;
  w = w.*self.bases{q}.weight(x(:,q));
end
