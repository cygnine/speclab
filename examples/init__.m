function[examples] = init__()
% init__ -- Initialization file for speclab/examples module
%
% [nodes] = init__()

module_list = {'chebyshev', 'legendre', 'jacobi', 'fourier', 'wiener'};

examples = recurse_files(pwd, module_list);

%examples.chebyshev = matlab_import('chebyshev');
%examples.legendre = matlab_import('legendre');
%examples.jacobi = matlab_import('jacobi');
%examples.fourier = matlab_import('fourier');
%examples.wiener = matlab_import('wiener');
