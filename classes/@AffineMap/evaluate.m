function[c] = evaluate(self, x)
% evaluate -- Evaluation method for AffineMap
%
% c = evaluate(self, x)
%
%     Evaluates the affine map at locations x in R^m. This function attempts to
%     determine which array dimension of the input x corresponds to the domain
%     dimension. In instances of conflicting possibilities, the first dimension
%     is prioritized.

xsize = size(x);
dim = find(xsize==self.m, 1, 'first');
if isempty(dim)
  error('Input array does not have compatible dimensions');
end
ysize = circshift(xsize, [0 1-dim]);
if length(ysize)>1
  ysize = ysize(2:end);
else
  ysize = 1;
end
y = [shiftdim(x, dim-1); ones([1 ysize])];

c = self.Aaug*y;
c = shiftdim(c(1:self.n,:), length(xsize) - (dim-1));
