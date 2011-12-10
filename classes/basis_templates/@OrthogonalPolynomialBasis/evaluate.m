function[p] = evaluate(self,x,n,varargin)
% p = evaluate(x,n,{d=0, normalization=self.normalization})
%
%     Evaluates the d-th derivative of the orthogonal polynomial basis at the
%     locations x. Each polynomial of linear index n is evaluated. p is a
%     length(x(:)) x length(n(:)) array. One may override the class instance's
%     normalization by the optional input 'normalization' (which must be of
%     FunctionNormalization derived type). 
%
%     If self.dim > 1, p is a size(x,1) x length(n(:)) array.

persistent input_parser 
if isempty(input_parser)
  from labtools import input_parser

end

params = {'d', 'normalization'};
defaults = {0, self.normalization};

opt = input_parser(params, defaults, [], varargin{:});
%inparse = inputParser();
%inparse.KeepUnmatched = true;
%
%inparse.addParamValue('d', 0);
%inparse.addParamValue('normalization', self.normalization);
%
%inparse.parse(varargin{:});
%opt = inparse.Results;

[n_array, nsize, numeln] = self.indexing(n);
N = max(n_array)+2;

xsize = size(x);  % reshape at end
x = x(:);

x = self.map_to_standard_domain(x.').';
[a,b] = self.recurrence(0:N);
p = self.evaluate_driver(x(:),a,b,n_array(:), opt.d(:));

% self.scale_functions also calls the recurrence formula -- any way to remove
% this waste of computation?
for q = 1:length(opt.d(:))
  p(:,:,q) = self.scale_functions(p(:,:,q),n_array,opt.normalization);
end

if size(p,3)==1
  if numeln == 1
    p = reshape(p, xsize);
  elseif numel(x) == 1
    p = reshape(p, nsize);
  end
end
