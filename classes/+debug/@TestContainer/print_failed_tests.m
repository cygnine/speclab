function[] = print_failed_tests(self)
% PRINT_FAILED_TESTS(SELF)
%
%     Pretty-prints information for each test that failed.

any_failures = false;

for n = 1:self.N
  test = self.tests{n};
  if ~test.result
    test.print()
    any_failures = true;
  end
end

if ~any_failures
  fprintf('All tests passed\n')
end
