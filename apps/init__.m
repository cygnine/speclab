function[applications] = init__()
% init__ -- Initialization file for speclab/applications module
%
% [nodes] = init__()

module_list = {'stein', 'genz_functions'};

applications = recurse_files(pwd, module_list);
%applications = add_module(applications, module_list);

%applications.stein = matlab_import('stein');
%applications.genz_functions = matlab_import('genz_functions');
