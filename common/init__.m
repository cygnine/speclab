function[common] = init__()
% init__ -- Initialization file for speclab/common module
%
% [nodes] = init__()

module_list = {'tensor'};
%common = recurse_files(pwd, module_list);

common.module_list = module_list;
common.recurse_files = true;
common.addpaths = {};
