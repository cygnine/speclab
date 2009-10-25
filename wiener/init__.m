function[wiener] = init__()
% init__ -- Initialization file for speclab/wiener module
%
% [nodes] = init__()

wiener = recurse_files(pwd);

wiener.quad = matlab_import('quad');
wiener.eval = matlab_import('eval');
fourier.maps = matlab_import('maps');
wiener.weights = matlab_import('weights');
wiener.coefficients = matlab_import('coefficients');
wiener.matrices = matlab_import('matrices');
wiener.fft = matlab_import('fft');
