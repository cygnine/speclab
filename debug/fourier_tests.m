function[container] = fourier_tests()
% [container] = fourier_tests()
%
%     Returns a TestContainer with Fourier tests.

fprintf('Building Fourier validation tests....\n');

global handles;
fourier = handles.speclab.fourier;
classic_tests = handles.speclab.debug.classic_fourier_tests;
approx_tests = handles.speclab.debug.fourier_approximation_tests;

import debug.*

container = TestContainer();

%%%% Regular Fourier case
opt = fourier.defaults();
opt.N = 10 + ceil(90*rand());
container = classic_tests(container,opt);
temp = opt.N;
opt.N = 80 + ceil(20*rand());
container = approx_tests(container,opt);
opt.N = temp;

opt.scale = 3*rand();
opt.shift = randn();
container = classic_tests(container,opt);
temp = opt.N;
opt.N = 80 + ceil(20*rand());
container = approx_tests(container,opt);
opt.N = temp;
