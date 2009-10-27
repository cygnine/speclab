function[tests] = all_tests()
% [TESTS] = ALL_TESTS()
%
%     Runs all validation tests for speclab.

debug = from_package_import_as('speclab', 'debug');
%debug = packages.speclab.debug;
tests = debug.linalg_tests();
tests = tests.extend(debug.jacobi_tests());
tests = tests.extend(debug.laguerre_tests());
tests = tests.extend(debug.fourier_tests());
tests = tests.extend(debug.wiener_tests());

tests = tests.run_tests();
tests.print_failed_tests('brief');
