function[orthopoly] = init__()
% init__ -- Initialization file for speclab/orthopoly module
%
% [nodes] = init__()

module_list = {'jacobi', 'interp', 'hermite', 'laguerre'};
%orthopoly = recurse_files(pwd, module_list);

orthopoly.module_list = module_list;
orthopoly.recurse_files = true;
orthopoly.addpaths = {};
