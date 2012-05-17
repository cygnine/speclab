function[theta,w] = stroud2(varargin)
% stroud2 -- Stroud's optimal degree-2 interpolation formula
%
% [theta,w] = stroud2([dim=1, scale=1, shift=0])
%
%     Computes the Stroud-2 quadrature formula in dim dimensions, which exactly
%     integrates polynomials of total degree <= 2. The choice of points is not
%     unique and is simply as regurgitation of the example given in [1], which
%     produces interior points on [-1,1]^dim (for default arguments).
%
%     theta is a (dim+1) x dim matrix, where each row is point in
%     dim-dimensional space. w is a (dim+1) x 1 vector that contains the weights
%     for each node.
%
%     scale and shift are length-dim vectors containing the affine scale and
%     shift for each dimension. (see speclab.common.box_to_scaleshift)
%
%  [1]: "Remarks on the disposition of points in numerical integration
%       formulas", A.H. Stroud

persistent input_parser affine_transformation parser
if isempty(input_parser)
  from labtools import input_parser
  from speclab.common import affine_transformation

  [opt, parser] = input_parser({'dim', 'scale', 'shift'}, ...
                               {1, 1, 0}, ...
                               [], ...
                               varargin{:});
else
  parser.parse(varargin{:});
  opt = parser.Results;
end

%opt = strict_inputs({'dim', 'scale', 'shift'}, {1, 1, 0}, [], varargin{:});

theta = zeros([opt.dim+1 opt.dim]);
w = zeros([opt.dim+1 1]);

ks = (0:opt.dim).';

for q = 1:opt.dim
  r = ceil(q/2);
  if mod(q,2)==1
    theta(:,q) = sqrt(2/3)*cos(2*r*ks*pi/(opt.dim+1));
  else
    theta(:,q) = sqrt(2/3)*sin(2*r*ks*pi/(opt.dim+1));
  end
end

if mod(opt.dim,2)==1
  theta(:,opt.dim) = sqrt(1/3)*(-1).^(ks);
end

w(:) = 2^opt.dim/(opt.dim+1);

theta = affine_transformation(theta, opt.scale, opt.shift);
w = w*prod(opt.scale);
