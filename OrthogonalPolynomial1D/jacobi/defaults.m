function[opt] = defaults(varargin);
% [OPT] = DEFAULTS(VARARGIN);
%     Function that generates default values for Jacobi polynomial expansions
%     given key-value pairs as input. The output is a struct of the results by
%     key-value combination.
%
% Defaults: alpha: -1/2                Jacobi polynomial parameter
%           beta: -1/2                 Jacobi polynomial parameter
%           scale: 1                   affine scale
%           shift: 0                   affine shift
%           normalization: 'normal'    Polynomial normalization
%           r: 1                       Radau point
%           r1: -1                     Lobatto point 1
%           r2: 1                      Lobatto point 2
% 
global handles;

jnames = {'alpha', 'beta', 'shift', 'scale','d','normalization','r','r1','r2'};
jdefaults = {-1/2, -1/2, 0, 1,0,'normal',1,-1,1};
opt = handles.common.InputSchema(jnames,jdefaults,[],varargin{:});
