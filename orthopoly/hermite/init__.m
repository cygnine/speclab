function[hermite] = init__()
% init__ -- Initialization file for speclab/orthopoly/hermite module
%
% [nodes] = init__()

hermite = recurse_files(pwd);

add_to_deprecation_list('speclab.orthopoly.hermite.eval', ...
                        'speclab.orthopoly.hermite.quad', ...
                        'speclab.orthopoly.hermite.coefficients', ...
                        'speclab.orthopoly.hermite.weight');
hermite.quad = matlab_import_deprecated('quad');
hermite.eval = matlab_import_deprecated('eval');
hermite.coefficients = matlab_import_deprecated('coefficients');
hermite.weights = matlab_import_deprecated('weights');
