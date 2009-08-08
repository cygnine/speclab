function[] = run_test(obj)
% VALIDATIONTEST.RUN_TEST()
%
%     Runs the ValidationTest and sticks the boolean result in OBJ.RESULT.

try
  data = obj.get_data;
  try
    obj.result = obj.validator(data);
  catch obj.error_exception
    obj.result = false;
  end
catch obj.error_exception
  obj.result = false;
end
