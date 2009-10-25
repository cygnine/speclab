classdef ValidationTest
  properties
    result = false; % Boolean result for test
    description = ''; % Short description of test
    validator = @() false;  % Validation function
    data_generator = @(x,opt) []; % Function generating data for validator
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

      global packages;
      inputs = {'description', 'validator', 'data_generator', 'parameters'};
      defaults = {'No description', @() false, @(x,opt) [], struct([])};
      opt = packages.labtools.input_schema(inputs, defaults, [], varargin{:});

      obj.description = opt.description;
      obj.parameters = opt.parameters;
      obj.validator = opt.validator;
      obj.data_generator = opt.data_generator;
    end

    data = get_data(obj); % Retrieves data created by data_generator
    obj = run_test(obj); % Runs the test and stores outcome in result
    print(obj,varargin); % Pretty-prints the test description and parameters
  end
end
