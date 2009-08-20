function[container] = wiener_tests()
% [container] = wiener_tests()
%
%     Returns a TestContainer with Wiener tests.

fprintf('Building Wiener validation tests....\n');

global handles;
wiener = handles.speclab.wiener;
classic_tests = handles.speclab.debug.classic_wiener_tests;
approx_tests = handles.speclab.debug.wiener_approximation_tests;
fft_tests = handles.speclab.debug.wiener_fft_tests;
matrix_tests = handles.speclab.debug.wiener_matrix_tests;

import debug.*
container = TestContainer();

%%%% Regular s=1 Wiener tests %%%%
opt = wiener.defaults();
opt.N = 10 + ceil(90*rand());
container = classic_tests(container,opt);
temp = opt.N;
opt.N = 80 + ceil(20*rand());
container = approx_tests(container,opt);
container = fft_tests(container,opt);
container = matrix_tests(container,opt);
opt.N = temp;

opt.scale = 3*rand();
opt.shift = randn();
container = classic_tests(container,opt);
temp = opt.N;
opt.N = 80 + ceil(20*rand());
container = approx_tests(container,opt);
container = fft_tests(container,opt);
container = matrix_tests(container,opt);
opt.N = temp;

%%%% Other s~=1 case (integer) %%%%
opt = wiener.defaults();
opt.N = 10 + ceil(90*rand());
opt.s = round(0.5 + 10*rand());
temp = opt.N;
opt.N = 80 + ceil(20*rand());
container = approx_tests(container,opt);
container = fft_tests(container,opt);
container = matrix_tests(container,opt);
opt.N = temp;

opt.scale = 3*rand();
opt.shift = randn();
temp = opt.N;
opt.N = 80 + ceil(20*rand());
container = approx_tests(container,opt);
container = fft_tests(container,opt);
container = matrix_tests(container,opt);
opt.N = temp;

%%%% Other s~=1 cases %%%%
opt = wiener.defaults();
opt.N = 10 + ceil(90*rand());
opt.s = 0.5 + 10*rand();
temp = opt.N;
opt.N = 80 + ceil(20*rand());
container = approx_tests(container,opt);
container = matrix_tests(container,opt);
opt.N = temp;

opt.scale = 3*rand();
opt.shift = randn();
temp = opt.N;
opt.N = 80 + ceil(20*rand());
container = approx_tests(container,opt);
container = matrix_tests(container,opt);
opt.N = temp;
