function[f] = continuous(x,varargin)
% continuous -- The "continuous" Genz test function
%
% f = continuous(x,{w=zeros([dim 1]),dim=size(x,2),c=zeros([dim 1])})
%
%     Evaluates the "continuous" Genz test function defined as
%
%        f(x) = \exp( - \sum_{i=1}^dim c_i |x_i - w_i| )
%
%     where x \in R^dim. x must be an N x dim matrix, where N is any number. If
%     dim is not specified, it is set to size(x,2). 

persistent strict_inputs
if isempty(strict_inputs)
  from labtools import strict_inputs
end

opt = strict_inputs({'dim', 'w', 'c'}, {size(x,2), [], []}, [], varargin{:});
if isempty(opt.c)
  opt.c = zeros([opt.dim 1]);
end
if isempty(opt.w)
  opt.w = zeros([opt.dim 1]);
end

f = abs(x - repmat(opt.w(:).', [size(x,1) 1]))*opt.c;
f = exp(-f);
