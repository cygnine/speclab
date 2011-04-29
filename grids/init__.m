function[grids] = init()
% init__ -- Initialization file for speclab.grids
%
% grids = init__()

module_list = {'qmc', 'tensor'};

grids = recurse_files(pwd, module_list);
%grids.qmc = matlab_import('qmc');
%grids.tensor = matlab_import('tensor');
