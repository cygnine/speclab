function[applications] = init__()
% init__ -- Initialization file for speclab/applications module
%
% [nodes] = init__()

module_list = {'stein', 'genz_functions'};
%applications = recurse_files(pwd, module_list);

applications.module_list = module_list;
applications.recurse_files = true;
applications.addpaths = {};
