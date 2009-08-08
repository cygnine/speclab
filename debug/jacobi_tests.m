function[container] = jacobi_tests()
% [CONTAINER] = JACOBI_TESTS()
%
%     Returns a TestContainer with all the Jacobi polynomial tests.

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
jquad_tests = handles.speclab.debug.jacobi_quadrature_tests_append;
import debug.*

container = TestContainer();

%%%%% Chebyshev case
opt = jac.defaults();
opt.N = 10 + ceil(90*rand());
container = jquad_tests(container,opt);

opt.shift = randn();
opt.scale = 3*rand();
container = jquad_tests(container,opt);
