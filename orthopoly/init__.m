function[orthopoly] = init__()
% init__ -- Initialization file for speclab/orthopoly module
%
% [nodes] = init__()

module_list = {'jacobi', 'interp', 'hermite', 'laguerre'};

orthopoly = recurse_files(pwd, module_list);

%orthopoly.jacobi = matlab_import('jacobi');
%orthopoly.interp = matlab_import('interp');
%orthopoly.hermite = matlab_import('hermite');
%orthopoly.laguerre = matlab_import('laguerre');
