function[] = print_failed_tests(self)
% PRINT_FAILED_TESTS(SELF)
%
%     Pretty-prints information for each test that failed.

for n = 1:self.N
  test = self.tests{n};
  if ~test.result
    test.print()
  end
end
