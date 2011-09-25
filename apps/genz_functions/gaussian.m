function[f] = gaussian(x,varargin)
% gaussian -- The "gaussian" Genz test function
%
% f = gaussian(x,{w=zeros([dim 1]),dim=size(x,2),c=ones([dim 1])})
%
%     Evaluates the "gaussian" Genz test function defined as
%
%        f(x) = \exp( - \sum_{i=1}^dim c_i^2 (x_i - w_i)^2 )
%
%     where x \in R^dim. x must be an N x dim matrix, where N is any number. If
%     dim is not specified, it is set to size(x,2). 

persistent strict_inputs
if isempty(strict_inputs)
  from labtools import strict_inputs
end

opt = strict_inputs({'dim', 'w', 'c'}, {size(x,2), [], []}, [], varargin{:});
if isempty(opt.c)
  opt.c = ones([opt.dim 1]);
end
if isempty(opt.w)
  opt.w = zeros([opt.dim 1]);
end
%disp([opt.w(:) opt.c(:)]);

f = (x - repmat(opt.w(:).', [size(x,1) 1])).^2*(opt.c(:).^2);
f = exp(-f);
