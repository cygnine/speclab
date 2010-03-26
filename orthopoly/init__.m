function[orthopoly] = init__()
% init__ -- Initialization file for speclab/orthopoly module
%
% [nodes] = init__()

orthopoly = recurse_files(pwd);

orthopoly.jacobi = matlab_import('jacobi');
orthopoly.interp = matlab_import('interp');
orthopoly.hermite = matlab_import('hermite');
orthopoly.laguerre = matlab_import('laguerre');
