function[interp] = init__()
% init__ -- Initialization file for speclab/orthopoly/interp module
%
% [nodes] = init__()

%interp = recurse_files(pwd);
%interp.least_utils = matlab_import('least_utils');

interp.module_list = {'least_utils'};
interp.recurse_files = true;
interp.addpaths = {};
