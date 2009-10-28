function[common] = init__()
% init__ -- Initialization file for speclab/common module
%
% [nodes] = init__()

common = recurse_files(pwd);
common.tensor = matlab_import('tensor');
