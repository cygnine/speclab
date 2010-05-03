function[f] = corner_peak(x,varargin)
% corner_peak -- The "corner peak" Genz test function
%
% f = corner_peak(x,{dim=size(x,2),c=ones([dim 1])})
%
%     Evaluates the "corner peak" Genz test function defined as
%
%        f(x) = (1 + \sum_{i=1}^dim c_i x_i )^{-d-1}
%
%     where x \in R^dim. x must be an N x dim matrix, where N is any number. If
%     dim is not specified, it is set to size(x,2). 

persistent strict_inputs
if isempty(strict_inputs)
  from labtools import strict_inputs
end

opt = strict_inputs({'dim', 'c'}, {size(x,2), []}, [], varargin{:});
if isempty(opt.c)
  opt.c = ones([opt.dim 1]);
end

f = (1 + x*opt.c).^(-opt.dim-1);
