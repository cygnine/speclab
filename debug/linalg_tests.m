function[container] = linalg_tests()
% [container] = linalg_tests()
%
%     Runs various tests for matlab-implemented linear algebra routines.

fprintf('Building linear algebra tests....\n');

import debug.*

%inversion_tests = packages.speclab.debug.linalg_inversion_tests;
from speclab.debug import linalg_inversion_tests as inversion_tests

container = TestContainer();

opt.N = 10 + ceil(40*rand());
opt.bandwidth = 1;
container = inversion_tests(container,opt);

opt.bandwidth = 2;
container = inversion_tests(container,opt);

opt.bandwidth = 4;
container = inversion_tests(container,opt);

opt.bandwidth = 7;
container = inversion_tests(container,opt);
