function[container] = jacobi_tests()
% [CONTAINER] = JACOBI_TESTS()
%
%     Returns a TestContainer with all the Jacobi polynomial tests.

fprintf('Building Jacobi polynomial validation tests....\n');

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
jquad_tests = handles.speclab.debug.jacobi_quadrature_tests_append;
japprox_tests = handles.speclab.debug.jacobi_approximation_tests_append;
jcoeff_tests = handles.speclab.debug.jacobi_coefficient_tests_append;
import debug.*

container = TestContainer();

%%%%% Chebyshev case
opt = jac.defaults();
opt.N = 10 + ceil(90*rand());
container = jquad_tests(container,opt);
container = japprox_tests(container,opt);
container = jcoeff_tests(container,opt);

opt.shift = randn();
opt.scale = 3*rand();
container = jquad_tests(container,opt);
container = japprox_tests(container,opt);
container = jcoeff_tests(container,opt);

%%%%% Legendre case
opt = jac.defaults();
opt.alpha = 0; opt.beta = 0;
opt.N = 10 + ceil(90*rand());
container = jquad_tests(container,opt);
container = japprox_tests(container,opt);
container = jcoeff_tests(container,opt);

opt.shift = randn();
opt.scale = 3*rand();
container = jquad_tests(container,opt);
container = japprox_tests(container,opt);
container = jcoeff_tests(container,opt);

%%%%% Random (alpha,beta) case
opt = jac.defaults();
opt.alpha = -1 + 10*rand(); opt.beta = -1 + 10*rand();
opt.N = 10 + ceil(90*rand());
container = jquad_tests(container,opt);
container = japprox_tests(container,opt);
container = jcoeff_tests(container,opt);

opt.shift = randn();
opt.scale = 3*rand();
container = jquad_tests(container,opt);
container = japprox_tests(container,opt);
container = jcoeff_tests(container,opt);

%%%%% Random (alpha,beta) case (yes, again)
opt = jac.defaults();
opt.alpha = -1 + 10*rand(); opt.beta = -1 + 10*rand();
opt.N = 10 + ceil(90*rand());
container = jquad_tests(container,opt);
container = japprox_tests(container,opt);
container = jcoeff_tests(container,opt);

opt.shift = randn();
opt.scale = 3*rand();
container = jquad_tests(container,opt);
container = japprox_tests(container,opt);
container = jcoeff_tests(container,opt);
