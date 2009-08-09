classdef TestContainer
  properties
    N = 0; % The number of tests in the container
    tests = cell(0); % Cell array of tests
    failed_test_indices = [];  % Linear indices of failed tests
    failed_tests = {};
  end
  methods

    function self = TestContainer()
    % OBJ = TESTCONTAINER()
    %
    %     Initializes an object that holds various ValidationTests. Use the
    %     .append or .extend methods to add ValidationTest(s).
      self.N = 0;
      tests = cell(0);
    end

    self = append(self,test); % Adds test (type ValidationTest) to the container
    self = extend(self,other); % Concatenates tests in other (type TestContainer) to self
    self = run_tests(self); % Runs all the test in the container
    function tests = get.failed_tests(self)
      tests = self.tests(self.failed_test_indices);
    end
    print_failed_tests(self,varargin); % Prints information for all failed tests
  end
end
