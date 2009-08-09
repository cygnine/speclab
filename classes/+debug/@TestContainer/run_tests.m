function self = run_tests(self)
% SELF = RUN_TESTS(SELF)
%
%     Runs all the ValidationTests contained in SELF. This does not print any
%     output; all the results are stored in each individual test object's data.

fprintf('     %d tests to run\n', self.N);
fprintf('%%%%%%%%       Legend       %%%%%%%%\n     * : Passed test\n     F : Failed test\n\n');
fprintf('       ');

for n = 1:self.N
  self.tests{n} = self.tests{n}.run_test();
  if self.tests{n}.result
    fprintf('*');
  else 
    fprintf('F');
    self.failed_test_indices = [self.failed_test_indices, n];
  end
  if mod(n,10)==0
    fprintf('\n      ')
    if mod(n,50)==0
      fprintf('------------ %d tests completed\n       ', n)
    else
      fprintf(' ');
    end
  end
end

fprintf('\n');
