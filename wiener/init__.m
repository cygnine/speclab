function[wiener] = init__()
% init__ -- Initialization file for speclab/wiener module
%
% [nodes] = init__()

module_list = {'quad', 'eval', 'maps', 'weights', 'coefficients', 'matrices', ...
               'operators', 'fft'};
%wiener = recurse_files(pwd, module_list);

wiener.module_list = module_list;
wiener.recurse_files = true;
wiener.addpaths = {};
