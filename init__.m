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
speclab.newton_polynomials = matlab_import('newton_polynomials');

speclab.debug = matlab_import('debug');
speclab.examples = matlab_import('examples');
speclab.apps = matlab_import('apps');
speclab.quad = matlab_import('quad');
speclab.filter = matlab_import('filter');

pwd_addpath('classes');

speclab.orthopoly1d = speclab.orthopoly;
add_to_deprecation_list('speclab.orthopoly1d');
