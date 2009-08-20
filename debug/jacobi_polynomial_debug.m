function[varargout] = jacobi_polynomial_debug()
% [FLAGS,DESCRIPTIONS,PARAMETERS] = JACOBI_POLYNOMIAL_DEBUG
% 
%     Runs various tests on speclab's Jacobi polynomial package. FLAGS is a
%     boolean array containing success indicators for each test, whose
%     description is in the cell array DESCRIPTIONS. PARAMETERS is a cell vector
%     with each element holding a struct with the parameters for the expansion.

fprintf('Running Jacobi polynomial tests...');

global handles;
jac = handles.speclab.orthopoly1d.jacobi;
debug = handles.speclab.debug;
dcat = debug.concatenate_debug_information;

flags = false(0);
descriptions = cell(0);
parameters = cell(0);

%%%%%%%%%%% Chebyshev case %%%%%%%%%%% 
N = 20 + ceil(90*rand);
opt = jac.defaults('N',N);
[temp,temp2] = debug.jacobi_polynomial_driver(opt);
temp2 = debug.description_prepend(temp2,'Chebyshev polynomial');
[flags,descriptions,parameters] = dcat(flags,descriptions,parameters, ...
   temp, temp2, opt);

%%%%%%%%%%% Chebyshev case, shifted, scaled %%%%%%%%%%% 
N = 20 + ceil(90*rand);
scale = 3*rand;
shift = 10*randn;
opt = jac.defaults('N',N,'scale',scale,'shift',shift);
[temp,temp2] = debug.jacobi_polynomial_driver(opt);
temp2 = debug.description_prepend(temp2,'Chebyshev polynomial, shift+scale');
[flags,descriptions,parameters] = dcat(flags,descriptions,parameters, ...
   temp, temp2, opt);

%%%%%%%%%%% Legendre case %%%%%%%%%%% 
N = 20 + ceil(90*rand);
opt = jac.defaults('N',N,'alpha',0,'beta',0);
[temp,temp2] = debug.jacobi_polynomial_driver(opt);
temp2 = debug.description_prepend(temp2,'Legendre polynomial');
[flags,descriptions,parameters] = dcat(flags,descriptions,parameters, ...
   temp, temp2, opt);

%%%%%%%%%%% Legendre case, shifted,scaled %%%%%%%%%%% 
N = 20 + ceil(90*rand);
scale = 3*rand;
shift = 10*randn;
opt = jac.defaults('N',N,'alpha',0,'beta',0,'scale',scale,'shift',shift);
[temp,temp2] = debug.jacobi_polynomial_driver(opt);
temp2 = debug.description_prepend(temp2,'Legendre polynomial, scale+shift');
[flags,descriptions,parameters] = dcat(flags,descriptions,parameters, ...
   temp, temp2, opt);

%%%%%%%%%%% Random (a,b) case %%%%%%%%%%% 
N = 20 + ceil(90*rand);
alpha = -1 + 10*rand;
beta = -1 + 10*rand;
opt = jac.defaults('N',N,'alpha',alpha,'beta',beta);
[temp,temp2] = debug.jacobi_polynomial_driver(opt);
temp2 = debug.description_prepend(temp2,'Asymmetric Jacobi polynomial');
[flags,descriptions,parameters] = dcat(flags,descriptions,parameters, ...
   temp, temp2, opt);

%%%%%%%%%%% Random (a,b) case, shift+scale %%%%%%%%%%% 
N = 20 + ceil(90*rand);
alpha = -1 + 10*rand;
beta = -1 + 10*rand;
scale = 3*rand;
shift = 10*randn;
opt = jac.defaults('N',N,'alpha',alpha,'beta',beta,'scale',scale,'shift',shift);
[temp,temp2] = debug.jacobi_polynomial_driver(opt);
temp2 = debug.description_prepend(temp2,'Asymmetric Jacobi polynomial, shift+scale');
[flags,descriptions,parameters] = dcat(flags,descriptions,parameters, ...
   temp, temp2, opt);

varargout{1} = flags;
varargout{2} = descriptions;
varargout{3} = parameters;

fprintf('done\n');
if any(not(flags))
  fprintf('    Some tests failed\n')
else
  fprintf('    All tests passed\n')
end
