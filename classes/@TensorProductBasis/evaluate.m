function[V] = evaluate(self,x,n, varargin)
% evaluate -- Evaluates tensor-product basis
%
% V = evaluate(self, x, n)
%
%     Evaluates the (one-based-indexed) n'th basis functions at the locations x.
%     x is an array of size N x self.dim, where each row is a point in
%     d-dimensional space. The output V is an array of size N x length(n), where
%     each column represents the appropriate basis function n evaluated at all
%     locations x.

persistent input_parser parser
if isempty(parser)
  from labtools import input_parser
  [opt, parser] = input_parser({'d', 'normalization'}, ...
                               {[], self.normalization}, ...
                               [], ...
                               varargin{:});
else
  parser.parse(varargin{:});
  opt = parser.Results;
end

if size(x,2) ~= self.dim
  error('Input points must have the same dimension (columns) as self.dim');
end
if not(isempty(opt.d)) && size(opt.d,2) ~= self.dim
  error('Specification of partial derivatives must have the same dimension (columns) as self.dim');
end

[n_array, nsize, numeln] = self.indexing(n);

if isempty(opt.d)
  Vsize = [size(x,1), numeln];
else
  Vsize = [size(x, 1), numeln, size(opt.d,1)];
end
V = ones(Vsize);

for q = 1:self.dim
  if isempty(opt.d)
    V = V.*reshape(self.bases{q}(x(:,q), n_array(:,q)), Vsize);
  else
    V = V.*reshape(self.bases{q}(x(:,q), n_array(:,q), 'd', opt.d(:,q)), Vsize);
  end
end

if (size(V,3)==1) && (size(x,1)==1)
  V = reshape(V, nsize);
end
