function[opt] = defaults(varargin);
% [opt] = defaults({s=0,t=0,scale=1,shift=0,x=false,n=false});
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

jnames = {'s', 't', 'shift', 'scale','x','n'};
jdefaults = {1, 0, 0, 1, false, false};
opt = handles.common.input_schema(jnames,jdefaults,[],varargin{:});
