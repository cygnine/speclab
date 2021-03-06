function[container] = laguerre_tests()
% [container] = laguerre_tests()
%
%     Returns a TestContainer with all the Laguerre polynomial tests.

fprintf('Building Laguerre polynomial validation tests....\n');

from speclab.orthopoly import laguerre as lag
from speclab.debug import laguerre_quadrature_tests_append as lquad_tests
from speclab.debug import laguerre_approximation_tests_append as lapprox_tests
from speclab.debug import laguerre_coefficient_tests_append as lcoeff_tests
%[lquad_tests, lapprox_tests, lcoeff_tests] = from_package_import_as(...
%  'speclab.debug', 'laguerre_quadrature_tests_append', ...
%                   'laguerre_approximation_tests_append', ...
%                   'laguerre_coefficient_tests_append');
%lag = packages.speclab.orthopoly.laguerre;
%lquad_tests = packages.speclab.debug.laguerre_quadrature_tests_append;
%lapprox_tests = packages.speclab.debug.laguerre_approximation_tests_append;
%lcoeff_tests = packages.speclab.debug.laguerre_coefficient_tests_append;
%jop_tests = packages.speclab.debug.laguerre_operator_tests_append;
import debug.*

container = TestContainer();

opt = lag.defaults();
opt.N = 10 + ceil(90*rand());
container = lquad_tests(container,opt);
temp = opt.N;
opt.N = 70 + ceil(30*rand());
container = lapprox_tests(container,opt);
container = lcoeff_tests(container,opt);
%container = jop_tests(container,opt);

opt.shift = randn();
opt.scale = 3*rand();
container = lquad_tests(container,opt);
temp = opt.N;
opt.N = 70 + ceil(30*rand());
container = lapprox_tests(container,opt);
container = lcoeff_tests(container,opt);
%container = jop_tests(container,opt);

opt.alpha = -1 + 10*rand();
opt.scale = 1;
opt.shift = 0;
container = lquad_tests(container,opt);
temp = opt.N;
opt.N = 70 + ceil(30*rand());
container = lapprox_tests(container,opt);
container = lcoeff_tests(container,opt);

opt.shift = randn();
opt.scale = 3*rand();
container = lquad_tests(container,opt);
temp = opt.N;
opt.N = 70 + ceil(30*rand());
container = lapprox_tests(container,opt);
container = lcoeff_tests(container,opt);
