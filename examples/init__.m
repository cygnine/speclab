function[examples] = init__()
% init__ -- Initialization file for speclab/examples module
%
% [nodes] = init__()

module_list = {'chebyshev', 'legendre', 'jacobi', 'fourier', 'wiener'};
%examples = recurse_files(pwd, module_list);

examples.module_list = module_list;
examples.recurse_files = true;
examples.addpaths = {};
