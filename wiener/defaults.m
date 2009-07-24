function[opt] = defaults(varargin);
% [OPT] = DEFAULTS({S=0,T=0,SCALE=1,SHIFT=0,X=false,N=false});
%     Function that generates default values for generalized Wiener expansions
%     given key-value pairs as input. The output is a struct of the results by
%     key-value combination.
%
% Defaults: s: 0                   Generalized Fourier parameter, {any # > -1/2}
%           t: 0                   Generalized Fourier parameter, {any # > -1/2}
%           scale: 1                   affine scale, {any positive #}
%           shift: 0                   affine shift, {any real #}
%           x: false                   Placeholder for optional input in vandermonde routines
%           n: false                   Placeholder for optional input in vandermonde routines
% 
global handles;

jnames = {'s', 't', 'shift', ...
          'scale','x','n'};
jdefaults = {0, 0, 0, 1, false, false};
opt = handles.common.InputSchema(jnames,jdefaults,[],varargin{:});
