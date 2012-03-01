function[speclab] = init__()
% init__ -- Initialization file for speclab package
%
% [speclab] = init__()

module_list = {'common', 'orthopoly', 'fourier', 'wiener', 'monomials', 'grids', ...
               'newton_polynomials', 'quad', 'debug', 'examples', 'apps', 'filter', ...
               'utils', 'd1_utils'};
%speclab = recurse_files(pwd, module_list);
%pwd_addpath('classes');
%pwd_addpath('classes/normalizations');
%pwd_addpath('classes/basis_templates');
%pwd_addpath('classes/indexing');
%pwd_addpath('classes/indexsets');
%pwd_addpath('classes/misc');
%speclab.orthopoly1d = speclab.orthopoly;
%add_to_deprecation_list('speclab.orthopoly1d');

speclab.module_list = module_list;
speclab.recurse_files = true;
speclab.addpaths = {'classes', 'classes/normalizations', 'classes/basis_templates', ...
                    'classes/indexing', 'classes/indexsets', 'classes/misc'};
