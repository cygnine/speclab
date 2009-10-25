function[speclab] = init__()
% init__ -- Initialization file for speclab package
%
% [nodes] = init__()

speclab = recurse_files(pwd);
speclab.common = matlab_import('common');
speclab.orthopoly1d = matlab_import('orthopoly1d');
speclab.fourier = matlab_import('fourier');
speclab.wiener = matlab_import('wiener');

speclab.monomials = matlab_import('monomials');
speclab.newton_polynomials = matlab_import('newton_polynomials');

speclab.debub = matlab_import('debug');
speclab.examples = matlab_import('examples');

pwd_addpath('classes');
