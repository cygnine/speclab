function[p] = evaluate(self,x,n,varargin)
% p = evaluate(x,n,{d=0, normalization=[]})
%
%     Evaluates the d-th derivative of the orthogonal polynomial basis at the
%     locations x. Each polynomial of linear index n is evaluated. p is a
%     length(x(:)) x length(n(:)) array. One may override the class instance's
%     normalization by the optional input 'normalization' (which must be of
%     FunctionNormalization derived type). 
%
%     If self.dim > 1, p is a size(x,1) x length(n(:)) array.

persistent indexing strict_inputs opoly_evaluate spdiag
if isempty(indexing)
  from speclab.common.tensor import linear_to_array_indexing as indexing
  from labtools import strict_inputs
  from speclab.d1_utils import opoly_evaluate
  from labtools import spdiag
end

opt = strict_inputs({'d', 'normalization'}, {0, self.normalization}, [], varargin{:});

%if self.dimension==1
n_array = n(:); % no need to call linear_to_array_indexing
N = max(n)+2;

xsize = size(x);  % reshape at end
x = x(:);

y = (2*x.^2 - 1);
even_n_flags = mod(n_array,2)==0;
odd_n_flags = mod(n_array,2)==1;

even_n = n_array(even_n_flags);
odd_n = n_array(odd_n_flags);

p = ones([size(x,1), size(n_array,1)]);
x = self.map_to_standard_domain(x.').';

%[a,b] = self.recurrence(0:N);
p(:,even_n_flags) = self.evenbasis(y, even_n/2)*spdiag(self.constant_cn(even_n/2));
p(:,odd_n_flags) = spdiag(x)*self.oddbasis(y, (odd_n-1)/2)*spdiag(self.constant_cn((odd_n-1)/2, self.lambda, self.mu+1));

%p = p.*opoly_evaluate(x(:),a,b,n_array(:), opt.d(:));

% self.scale_functions also calls the recurrence formula -- any way to remove
% this waste of computation?
p = self.scale_functions(p,n,opt.normalization);

if numel(n) == 1
  if self.dimension == 1
    p = reshape(p, xsize);
  else
    p = reshape(p, ysize);
  end
end
