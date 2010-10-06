function[opt] = defaults(varargin);
% [OPT] = DEFAULTS(VARARGIN);
%     Function that generates default values for Jacobi polynomial expansions
%     given key-value pairs as input. The output is a struct of the results by
%     key-value combination.
%
% Defaults: alpha: -1/2                Jacobi polynomial parameter, {any # > -1}
%           beta: -1/2                 Jacobi polynomial parameter, {any # > -1}
%           scale: 1                   affine scale, {any positive #}
%           shift: 0                   affine shift, {any real #}
%           d: 0                       Derivative indicator, [0,1,2,...]
%           normalization: 'normal'    Polynomial normalization, {'normal', 'monic'}
%           r: 1                       Radau point, [-1,1]
%           r1: -1                     Lobatto point 1, [-1,1]
%           r2: 1                      Lobatto point 2, [-1,1]
%           x: false                   Placeholder for optional input in vandermonde routines
%           n: false                   Placeholder for optional input in vandermonde routines
%           dim: 1                     Spatial dimension of the tensor-product evaluation
%           weight_normalization: ''   Specifies the normalization of the weight function
%           associated_index: 0        Specifies the associated polynomial index
% 

persistent input_schema
if isempty(input_schema)
  from labtools import input_schema
end

jnames = {'alpha', 'beta', 'shift', ...
          'scale','d','normalization','r','r1','r2','x','n','dim','weight_normalization','associated_index'};
jdefaults = {-1/2, -1/2, 0, 1,0,'normal',1,-1,1,false, false, 1, '', 0};
opt = input_schema(jnames,jdefaults,[],varargin{:});

% Change r, r1, r2 to match scale+shift
opt.r2 = opt.scale+opt.shift;
opt.r = opt.r;
opt.r1 = -opt.scale+opt.shift;

% Change inputs to have correct dimensionality
if opt.dim>1
  if length(opt.scale)==1
    opt.scale = opt.scale*ones([opt.dim 1]);
    opt.shift = opt.shift*ones([opt.dim 1]);
  end
  if length(opt.alpha)==1
    opt.alpha = opt.alpha*ones([opt.dim 1]);
    opt.beta = opt.beta*ones([opt.dim 1]);
  end
end
