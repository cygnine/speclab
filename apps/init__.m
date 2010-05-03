function[applications] = init__()
% init__ -- Initialization file for speclab/applications module
%
% [nodes] = init__()

applications = recurse_files(pwd);
applications.stein = matlab_import('stein');
applications.genz_functions = matlab_import('genz_functions');
