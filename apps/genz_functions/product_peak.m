function[f] = product_peak(x,varargin)
% product_peak -- The "product peak" Genz test function
%
% f = product_peak(x,{w=zeros([dim 1]),dim=size(x,2),c=ones([dim 1])})
%
%     Evaluates the "product peak" Genz test function defined as
%
%        f(x) = \prod_{i=1}^dim (1/c^2_i + (x_i- w_i)^2)^{-1}
%
%     where x \in R^dim. x must be an N x dim matrix, where N is any number. If
%     dim is not specified, it is set to size(x,2). The constants c_i are
%     usually restricted in the values (0, \infty), but no checking is performing in
%     this code to ensure that.

persistent input_parser parser
if isempty(input_parser)
  from labtools import input_parser
end
[opt, parser] = input_parser({'dim', 'w', 'c'}, ...
                             {size(x,2), [], []}, ...
                             [], ...
                             varargin{:});
parser.parse(varargin{:});
opt = parser.Results;

%opt = strict_inputs({'dim', 'w', 'c'}, {size(x,2), [], []}, [], varargin{:});
if isempty(opt.c)
  opt.c = ones([opt.dim 1]);
end
if isempty(opt.w)
  opt.w = zeros([opt.dim 1]);
end

f = (x - repmat(opt.w(:).', [size(x,1) 1])).^2 + 1./(repmat(opt.c(:).', [size(x,1) 1])).^2;
f = prod(1./f, 2);
