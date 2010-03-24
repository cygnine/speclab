function[interp] = init__()
% init__ -- Initialization file for speclab/orthopoly1d/interp module
%
% [nodes] = init__()

interp = recurse_files(pwd);

interp.least_utils = matlab_import('least_utils');
