function[hermite] = init__()
% init__ -- Initialization file for speclab/orthopoly1d/hermite module
%
% [nodes] = init__()

hermite = recurse_files(pwd);

hermite.quad = matlab_import('quad');
hermite.eval = matlab_import('eval');
hermite.coefficients = matlab_import('coefficients');
hermite.weights = matlab_import('weights');
