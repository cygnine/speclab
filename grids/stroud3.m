function[theta,w] = stroud3(varargin)
% stroud3 -- Stroud's optimal degree-3 interpolation formula
%
% [theta,w] = stroud2([dim=1, scale=1, shift=0])
%
%     Computes the Stroud-3 quadrature formula in dim dimensions, which exactly
%     integrates polynomials of total degree <= 3. The choice of points is not
%     unique and is simply as regurgitation of the example given in [1], which
%     produces interior points on [-1,1]^dim (for default arguments).
%
%     theta is a 2*dim x dim matrix, where each row is point in
%     dim-dimensional space. w is a 2*dim x 1 vector that contains the weights
%     for each node.
%
%     scale and shift are length-dim vectors containing the affine scale and
%     shift for each dimension. (see speclab.common.box_to_scaleshift)
%
%  [1]: "Remarks on the disposition of points in numerical integration
%       formulas", A.H. Stroud

persistent strict_inputs affine_transformation
if isempty(strict_inputs)
  from labtools import strict_inputs
  from speclab.common import affine_transformation
end

opt = strict_inputs({'dim', 'scale', 'shift'}, {1, 1, 0}, [], varargin{:});

theta = zeros([2*opt.dim opt.dim]);
w = zeros([2*opt.dim 1]);

ks = (1:2*opt.dim).';

for q = 1:opt.dim
  r = ceil(q/2);
  if mod(q,2)==1
    theta(:,q) = sqrt(2/3)*cos((2*r-1)*ks*pi/opt.dim);
  else
    theta(:,q) = sqrt(2/3)*sin((2*r-1)*ks*pi/opt.dim);
  end
end

if mod(opt.dim,2)==1
  theta(:,opt.dim) = sqrt(1/3)*(-1).^(ks);
end

w(:) = 2^opt.dim/(2*opt.dim);

theta = affine_transformation(theta, opt.scale, opt.shift);
w = w*prod(opt.scale);
