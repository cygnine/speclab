function[] = run_tests(self)
% RUN_TESTS(SELF)
%
%     Runs all the ValidationTests contained in SELF. This does not print any
%     output; all the results are stored in each individual test object's data.

for n = 1:self.N
  test = self.tests{n};
  test.run_test();
  self.tests{n} = test;
end
