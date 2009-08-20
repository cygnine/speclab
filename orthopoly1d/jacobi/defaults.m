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
%           normalization: 'normal'    Polynomial normalization, {'normal', 'monic'}
%           d: 0                       Derivative indicator, [0,1,2,...]
%           r: 1                       Radau point, [-1,1]
%           r1: -1                     Lobatto point 1, [-1,1]
%           r2: 1                      Lobatto point 2, [-1,1]
%           x: false                   Placeholder for optional input in vandermonde routines
%           n: false                   Placeholder for optional input in vandermonde routines
% 
global handles;

jnames = {'alpha', 'beta', 'shift', ...
          'scale','d','normalization','r','r1','r2','x','n'};
jdefaults = {-1/2, -1/2, 0, 1,0,'normal',1,-1,1,false, false};
opt = handles.common.InputSchema(jnames,jdefaults,[],varargin{:});
