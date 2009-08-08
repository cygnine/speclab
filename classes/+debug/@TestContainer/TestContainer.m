classdef TestContainer
  properties
    N = 0; % The number of tests in the container
    tests = cell(0); % Cell array of tests
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

    append(self,test); % Adds test (type ValidationTest) to the container
    extend(self,other); % Concatenates tests in other (type TestContainer) to self
    run_tests(self); % Runs all the test in the container
    print_failed_tests(self); % Prints information for all failed tests
  end
end
