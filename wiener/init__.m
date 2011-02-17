function[wiener] = init__()
% init__ -- Initialization file for speclab/wiener module
%
% [nodes] = init__()

module_list = {'quad', 'eval', 'maps', 'weights', 'coefficients', 'matrices', ...
               'operators', 'fft'};

wiener = recurse_files(pwd, module_list);

%wiener.quad = matlab_import('quad');
%wiener.eval = matlab_import('eval');
%wiener.maps = matlab_import('maps');
%wiener.weights = matlab_import('weights');
%wiener.coefficients = matlab_import('coefficients');
%wiener.matrices = matlab_import('matrices');
%wiener.operators = matlab_import('operators');
%wiener.fft = matlab_import('fft');
