function[laguerre] = init__()
% init__ -- Initialization file for speclab/orthopoly/laguerre module
%
% [nodes] = init__()

laguerre = recurse_files(pwd);

add_to_deprecation_list('speclab.orthopoly.laguerre.eval', ...
                        'speclab.orthopoly.laguerre.quad', ...
                        'speclab.orthopoly.laguerre.coefficients', ...
                        'speclab.orthopoly.laguerre.weights', ...
                        'speclab.orthopoly.laguerre.connection', ...
                        'speclab.orthopoly.laguerre.operators');
laguerre.quad = matlab_import_deprecated('quad');
laguerre.eval = matlab_import_deprecated('eval');
laguerre.coefficients = matlab_import_deprecated('coefficients');
laguerre.weights = matlab_import_deprecated('weights');
laguerre.connection = matlab_import_deprecated('connection');
laguerre.operators = matlab_import_deprecated('operators');
