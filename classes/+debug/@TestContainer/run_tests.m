function self = run_tests(self)
% SELF = RUN_TESTS(SELF)
%
%     Runs all the ValidationTests contained in SELF. This does not print any
%     output; all the results are stored in each individual test object's data.

for n = 1:self.N
  self.tests{n} = self.tests{n}.run_test();
end
