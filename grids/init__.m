function[grids] = init()
% init__ -- Initialization file for speclab.grids
%
% grids = init__()

module_list = {'qmc', 'tensor'};
%grids = recurse_files(pwd, module_list);

grids.module_list = module_list;
grids.recurse_files = true;
grids.addpaths = {};
