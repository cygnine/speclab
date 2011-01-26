function[speclab] = init__()
% init__ -- Initialization file for speclab package
%
% [nodes] = init__()

speclab = recurse_files(pwd);
speclab.common = matlab_import('common');
speclab.orthopoly = matlab_import('orthopoly');
speclab.fourier = matlab_import('fourier');
speclab.wiener = matlab_import('wiener');

speclab.monomials = matlab_import('monomials');
speclab.grids = matlab_import('grids');
speclab.newton_polynomials = matlab_import('newton_polynomials');

speclab.debug = matlab_import('debug');
speclab.examples = matlab_import('examples');
speclab.apps = matlab_import('apps');
speclab.filter = matlab_import('filter');
speclab.utils = matlab_import('utils');
speclab.d1_utils = matlab_import('d1_utils');

pwd_addpath('classes');

speclab.orthopoly1d = speclab.orthopoly;
add_to_deprecation_list('speclab.orthopoly1d');
