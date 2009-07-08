% Script that defaults generates default values for Jacobi polynomial
% expansions. The defaults are arranged in cells in key-value-like pairs for
% input into InputSchema.
% Interval is [shift-scale, shift+scale]

% Defaults: alpha = -1/2, beta=-1/2, scale=1, shift=0, 
%           normalization=normal, r=+1,
%           r1=-1, r2=+1
jnames = {'alpha', 'beta', 'shift', 'scale','d','normalization','r','r1','r2'};
jdefaults = {-1/2, -1/2, 0, 1,0,'normal',1,-1,1);
opt = handles.common.InputSchema(jnames,jdefaults,[],varargin{:});
