function[common] = init__()
% init__ -- Initialization file for speclab/common module
%
% [nodes] = init__()

module_list = {'tensor'};

common = recurse_files(pwd, module_list);
%common = add_module(common, module_list);
%common.tensor = matlab_import('tensor');
