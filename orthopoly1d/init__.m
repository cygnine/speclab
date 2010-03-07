function[orthopoly1d] = init__()
% init__ -- Initialization file for speclab/orthopoly1d module
%
% [nodes] = init__()

orthopoly1d = recurse_files(pwd);

orthopoly1d.jacobi = matlab_import('jacobi');
orthopoly1d.interp = matlab_import('interp');
orthopoly1d.hermite = matlab_import('hermite');
orthopoly1d.laguerre = matlab_import('laguerre');
