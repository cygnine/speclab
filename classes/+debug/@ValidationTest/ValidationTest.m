classdef ValidationTest
  properties
    result = false; % Boolean result for test
    description = ''; % Short description of test
    validator = @() false;  % Validation function
    data_generator = @(x) []; % Function generating data for validator
    parameters = struct([]); % Parameters for generating data
    error_exception = MException('ValidationTest:OK', 'No error');
  end
  methods

    function obj = ValidationTest(varargin)
    % OBJ = VALIDATIONTEST({DESCRIPTION='NO DESCRIPTION', VALIDATOR=@() FALSE, 
    %                       DATA_GENERATOR=@(X) [], PARAMETERS=STRUCT([])})
    %
    %     Constructor method for classtype ValidationTest. This object is a
    %     boolean test to validate other code. 

      global handles;
      inputs = {'description', 'validator', 'data_generator', 'parameters'};
      defaults = {'No description', @() false, @(x) [], struct([])};
      opt = handles.common.InputSchema(inputs, defaults, [], varargin{:});

      obj.description = opt.description;
      obj.parameters = opt.parameters;
      obj.validator = opt.validator;
      obj.data_generator = opt.data_generator;
    end

    data = get_data(obj); % Retrieves data created by data_generator
    run_test(obj); % Runs the test and stores outcome in result
    print(obj); % Pretty-prints the test description and parameters
  end
end
