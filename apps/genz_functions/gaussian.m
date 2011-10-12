function[f] = gaussian(x,varargin)
% gaussian -- The "gaussian" Genz test function
%
% f = gaussian(x,{w=zeros([dim 1]),dim=size(x,2),c=ones([dim 1]), d=0})
%
%     Evaluates the "gaussian" Genz test function defined as
%
%        f(x) = \exp( - \sum_{i=1}^dim c_i^2 (x_i - w_i)^2 )
%
%     where x \in R^dim. x must be an N x dim matrix, where N is any number. If
%     dim is not specified, it is set to size(x,2). 
%
%     If the optional input d is equal to 0, then the function is evaluated. If
%     d is a sized (dim x 1) vector, then the directional derivative in the
%     direction d is computed and returned.

persistent input_parser parser
if isempty(parser)
  from labtools import input_parser
  [opt,parser] = input_parser({'dim', 'w', 'c', 'd'}, ...
                              {[], [], [], 0}, ...
                              [], ...
                              varargin{:});
else
  parser.parse(varargin{:});
  opt = parser.Results;
end

if isempty(opt.dim)
  opt.dim = size(x,2);
end
if isempty(opt.c)
  opt.c = ones([opt.dim 1]);
end
if isempty(opt.w)
  opt.w = zeros([opt.dim 1]);
end
%disp([opt.w(:) opt.c(:)]);

fexp = (x - repmat(opt.w(:).', [size(x,1) 1])).^2*(opt.c(:).^2);
f = exp(-fexp);
if isempty(opt.d) || all(opt.d == 0)
  return;
elseif size(opt.d)==[1 opt.dim];
  % Compute directional derivative
  % Make sure d is normalized
  opt.d = opt.d/norm(opt.d);

  linear_factors = 2*(x - repmat(opt.w(:).', [size(x,1) 1]))*diag(opt.c(:).^2);
  f = -exp(-f).*(linear_factors*opt.d(:));
else
  error('Unrecognized shape or type for optional input ''d''');
end
