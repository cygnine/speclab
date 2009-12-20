function[opt] = defaults(varargin);
% [opt] = defaults(varargin);
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

persistent input_schema
if isempty(input_schema)
  from labtools import input_schema
end

jnames = {'gamma', 'delta', 'shift', ...
          'scale','normalization','x','n'};
jdefaults = {0, 0, 0, 1,0,'normal',false, false};
opt = input_schema(jnames,jdefaults,[],varargin{:});
