function[container] = linalg_tests()
% [container] = linalg_tests()
%
%     Runs various tests for matlab-implemented linear algebra routines.

fprintf('Building linear algebra tests....\n');

global handles;
import debug.*

inversion_tests = handles.speclab.debug.linalg_inversion_tests;

container = TestContainer();

opt.N = 10 + ceil(40*rand());
opt.bandwidth = 1;
container = inversion_tests(container,opt);

opt.bandwidth = 2;
container = inversion_tests(container,opt);

opt.bandwidth = 5;
container = inversion_tests(container,opt);

opt.bandwidth = 10;
container = inversion_tests(container,opt);
