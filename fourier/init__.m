function[fourier] = init__()
% init__ -- Initialization file for speclab/fourier module
%
% [nodes] = init__()

fourier = recurse_files(pwd);

fourier.quad = matlab_import('quad');
fourier.eval = matlab_import('eval');
fourier.maps = matlab_import('maps');
fourier.weights = matlab_import('weights');
fourier.connection = matlab_import('connection');
fourier.fft = matlab_import('fft');
