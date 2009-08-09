function[] = print_failed_tests(self,varargin)
% print_failed_tests(self,{fmt=[]})
%
%     Pretty-prints information for each test that failed.

any_failures = false;
fmt = varargin{1};

for n = 1:self.N
  test = self.tests{n};
  if ~test.result
    test.print(fmt)
    any_failures = true;
  end
end

if ~any_failures
  fprintf('All tests passed\n')
end
