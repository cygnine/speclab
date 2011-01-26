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

persistent indexing strict_inputs opoly_evaluate
if isempty(indexing)
  from speclab.common.tensor import linear_to_array_indexing as indexing
  from labtools import strict_inputs
  from speclab.d1_utils import opoly_evaluate
end

opt = strict_inputs({'d', 'normalization'}, {0, self.normalization}, [], varargin{:});

%if self.dimension==1
n_array = n(:); % no need to call linear_to_array_indexing
N = max(n)+2;

xsize = size(x);  % reshape at end
x = x(:);
%else
%  warning('This will soon be deprecated');
%  n_array = indexing(n,'dim',self.dimension);
%  N = max(max(n_array)) + 2;
%
%  % *sigh*, bean-counting to make two-dimensional array
%  % dimension 1: separate points
%  % dimension 2: domain dimension
%  xsize = size(x);
%  ind = find(xsize==self.dimension, 1, 'first');
%  x = shiftdim(x, ind-1);
%  ysize = circshift(xsize, [0, 1-ind]);
%  n2 = ysize(2);
%  %x = reshape(x, [self.dimension, n2]);
%  ysize = xsize; ysize(ind) = 1;
%end


p = ones([size(x,1), size(n_array,1)]);
x = self.map_to_standard_domain(x.').';
%for q = 1:self.dimension
[a,b] = self.recurrence(0:N);
p = p.*opoly_evaluate(x(:),a,b,n_array(:), opt.d(:));
  %p = p.*self.eval_driver(x(:,q),a,b,n_array(:,q), opt.d(:,q));
%end

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
