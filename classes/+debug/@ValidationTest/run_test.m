function self = run_test(self)
% SELF = VALIDATIONTEST.RUN_TEST(SELF)
%
%     Runs the ValidationTest and sticks the boolean result in SELF.RESULT.

try
  data = self.get_data;
  try
    self.result = self.validator(data,self.parameters);
  catch self.error_exception;
    self.result = false;
  end
catch self.error_exception;
  self.result = false;
end
