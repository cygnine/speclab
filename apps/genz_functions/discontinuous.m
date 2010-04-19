function[f] = discontinuous(x,varargin)
% discontinuous -- The "discontinuous" Genz test function
%
% f = discontinuous(x,{w=zeros([dim 1]),dim=size(x,2),c=zeros([dim 1])})
%
%     Evaluates the "discontinuous" Genz test function defined as
%
%        f(x) = \exp( \sum_{i=1}^dim c_i x_i ), assuming x_i \leq w_i
%             = 0,                              else
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

f = exp(x*opt.c);

flags = x > repmat(opt.w(:).', [size(x,1) 1]);
f(flags) = 0;
