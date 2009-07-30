function[varargout] = jacobi_polynomial_debug()
% [FLAGS,DESCRIPTIONS,PARAMETERS] = JACOBI_POLYNOMIAL_DEBUG
% 
%     Runs various tests on speclab's Jacobi polynomial package. FLAGS is a
%     boolean array containing success indicators for each test, whose
%     description is in the cell array DESCRIPTIONS. PARAMETERS is a cell vector
%     with each element holding a struct with the parameters for the expansion.

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
debug = handles.speclab.debug;
dcat = debug.concatenate_debug_information;

flags = [];
descriptions = cell(0);
parameters = cell(0);

%%%%%%%%%%% Chebyshev case %%%%%%%%%%% 
N = 10 + ceil(150*rand);
opt = jac.defaults('N',N);
[temp,temp2] = debug.jacobi_polynomial_driver(opt);
descriptions = debug.description_prepend(descriptions,'Chebyshev polynomial');
[flags,descriptions,parameters] = dcat(flags,descriptions,parameters, ...
   opt, temp, temp2);

%%%%%%%%%%% Chebyshev case, shifted, scaled %%%%%%%%%%% 
N = 10 + ceil(150*rand);
scale = 10*rand;
shift = 10*randn;
opt = jac.defaults('N',N,'scale',scale,'shift',shift);
[temp,temp2] = debug.jacobi_polynomial_driver(opt);
descriptions = debug.description_prepend(descriptions,'Chebyshev polynomial, shift+scale');
[flags,descriptions,parameters] = dcat(flags,descriptions,parameters, ...
   opt, temp, temp2);

%%%%%%%%%%% Legendre case %%%%%%%%%%% 
N = 10 + ceil(150*rand);
opt = jac.defaults('N',N,'alpha',0,'beta',0);
[temp,temp2] = debug.jacobi_polynomial_driver(opt);
descriptions = debug.description_prepend(descriptions,'Legendre polynomial');
[flags,descriptions,parameters] = dcat(flags,descriptions,parameters, ...
   opt, temp, temp2);

%%%%%%%%%%% Legendre case, shifted,scaled %%%%%%%%%%% 
N = 10 + ceil(150*rand);
scale = 10*rand;
shift = 10*randn;
opt = jac.defaults('N',N,'alpha',0,'beta',0,'scale',scale,'shift',shift);
[temp,temp2] = debug.jacobi_polynomial_driver(opt);
descriptions = debug.description_prepend(descriptions,'Legendre polynomial, scale+shift');
[flags,descriptions,parameters] = dcat(flags,descriptions,parameters, ...
   opt, temp, temp2);

%%%%%%%%%%% Random (a,b) case %%%%%%%%%%% 
N = 10 + ceil(150*rand);
alpha = -1 + 10*rand;
beta = -1 + 10*rand;
opt = jac.defaults('N',N,'alpha',alpha,'beta',beta);
[temp,temp2] = debug.jacobi_polynomial_driver(opt);
descriptions = debug.description_prepend(descriptions,'Asymmetric Jacobi polynomial');
[flags,descriptions,parameters] = dcat(flags,descriptions,parameters, ...
   opt, temp, temp2);

%%%%%%%%%%% Random (a,b) case, shift+scale %%%%%%%%%%% 
N = 10 + ceil(150*rand);
alpha = -1 + 10*rand;
beta = -1 + 10*rand;
scale = 10*rand;
shift = 10*randn;
opt = jac.defaults('N',N,'alpha',alpha,'beta',beta,'scale',scale,'shift',shift);
[temp,temp2] = debug.jacobi_polynomial_driver(opt);
descriptions = debug.description_prepend(descriptions,'Asymmetric Jacobi polynomial, shift+scale');
[flags,descriptions,parameters] = dcat(flags,descriptions,parameters, ...
   opt, temp, temp2);
