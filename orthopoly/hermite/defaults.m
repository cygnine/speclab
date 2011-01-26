function[opt] = defaults(varargin)
% defaults -- returns default parameters for Hermite polynomial expansions
%
% opt = defaults({mu=0, shift=0, scale=1, d=0, x=0, normalization='nomal'})
%     Returns the default parameters for Hermite polynomial expansions as fields
%     of a struct opt. Overwrites defaults values with any user-specified
%     inputs.
%
%     mu                  The parameter specifying the family of Hermite
%                         polynomials. Must be greater than or equal to 0.
%     shift               The affine shift (translation).
%     scale               The affine scale (dilation).
%     d                   An integer (or array of integers) greater than or
%                         equal to zero specifying which derivatives are to be
%                         evaluated.
%     x                   The location of the Gauss-Radau point.
%     normalization       'normal' -- the L^2 normalized polynomials
%                         'monic' -- the monic polynomials
%     dim: 1              Spatial dimension of the tensor-product evaluation
%     associated_index: 0 Specifies the associated polynomial index

persistent input_schema
if isempty(input_schema)
  from labtools import input_schema
end

hnames = {'mu', 'shift', 'scale', 'd', 'x','normalization','dim','associated_index',... 
          'weight_normalization'};
hdefaults = {0, 0, 1, 0, 0, 'normal',1, 0, ''};

opt = input_schema(hnames, hdefaults, [], varargin{:});

% Change default x to match shift, scale
opt.x = opt.x.*opt.scale + opt.shift;

% Change inputs to have correct dimensionality
if opt.dim>1
  if length(opt.scale)==1
    opt.scale = opt.scale*ones([opt.dim 1]);
    opt.shift = opt.shift*ones([opt.dim 1]);
  end
  if length(opt.mu)==1
    opt.mu = opt.mu*ones([opt.dim 1]);
  end
end
