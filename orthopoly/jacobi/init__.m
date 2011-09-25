function[jacobi] = init__()
% init__ -- Initialization file for speclab/orthopoly/jacobi module
%
% [nodes] = init__()

%jacobi = recurse_files(pwd);
%jacobi.fft = matlab_import('fft');

%add_to_deprecation_list('speclab.orthopoly.jacobi.eval', ...
%                        'speclab.orthopoly.jacobi.quad', ...
%                        'speclab.orthopoly.jacobi.coefficients', ...
%                        'speclab.orthopoly.jacobi.connection', ...
%                        'speclab.orthopoly.jacobi.weights', ...
%                        'speclab.orthopoly.jacobi.operators');

%jacobi.quad = matlab_import_deprecated('quad');
%jacobi.eval = matlab_import_deprecated('eval');
%jacobi.coefficients = matlab_import_deprecated('coefficients');
%jacobi.connection = matlab_import_deprecated('connection');
%jacobi.weights = matlab_import_deprecated('weights');
%jacobi.operators = matlab_import_deprecated('operators');

module_list = {'quad', 'eval', 'coefficients', 'connection', 'weights', 'operators', 'fft'};

jacobi.module_list = module_list;
jacobi.recurse_files = true;
jacobi.addpaths = {};
