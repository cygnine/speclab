function[jacobi] = init__()
% init__ -- Initialization file for speclab/orthopoly1d/jacobi module
%
% [nodes] = init__()

jacobi = recurse_files(pwd);

jacobi.quad = matlab_import('quad');
jacobi.eval = matlab_import('eval');
jacobi.coefficients = matlab_import('coefficients');
jacobi.connection = matlab_import('connection');
jacobi.weights = matlab_import('weights');
jacobi.operators = matlab_import('operators');
jacobi.fft = matlab_import('fft');
