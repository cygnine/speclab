function[f] = oscillatory(x,varargin)
% oscillatory -- The "oscillatory" Genz test function
%
% f = oscillatory(x,{w=0,dim=size(x,2),c=ones([dim 1])})
%
%     Evaluates the "oscillatory" Genz test function defined as
%
%        f(x) = \cos{ 2\pi w + \sum_{i=1}^dim c_i x_i }
%
%     where x \in R^dim. x must be an N x dim matrix, where N is any number. If
%     dim is not specified, it is set to size(x,2).

persistent strict_inputs
if isempty(strict_inputs)
  from labtools import strict_inputs
end

opt = strict_inputs({'dim', 'w', 'c'}, {size(x,2), 0, []}, [], varargin{:});
if isempty(opt.c)
  opt.c = ones([opt.dim 1]);
end

f = cos(2*pi*opt.w + x*opt.c);
