function[opt] = defaults(varargin);
% [OPT] = DEFAULTS(VARARGIN);
%     Function that generates default values for generalized Fourier expansions
%     given key-value pairs as input. The output is a struct of the results by
%     key-value combination.
%
% Defaults: gamma: 0                   Generalized Fourier parameter, {any # > -1/2}
%           delta: 0                   Generalized Fourier parameter, {any # > -1/2}
%           scale: 1                   affine scale, {any positive #}
%           shift: 0                   affine shift, {any real #}
%           normalization: 'normal'    Function normalization, {'normal'}
%           x: false                   Placeholder for optional input in vandermonde routines
%           n: false                   Placeholder for optional input in vandermonde routines
% 
global handles;

jnames = {'gamma', 'delta', 'shift', ...
          'scale','normalization','x','n'};
jdefaults = {0, 0, 0, 1,0,'normal',false, false};
opt = handles.common.InputSchema(jnames,jdefaults,[],varargin{:});
