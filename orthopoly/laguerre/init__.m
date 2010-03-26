function[laguerre] = init__()
% init__ -- Initialization file for speclab/orthopoly/laguerre module
%
% [nodes] = init__()

laguerre = recurse_files(pwd);

laguerre.quad = matlab_import('quad');
laguerre.eval = matlab_import('eval');
laguerre.coefficients = matlab_import('coefficients');
laguerre.weights = matlab_import('weights');
laguerre.connection = matlab_import('connection');
laguerre.operators = matlab_import('operators');
