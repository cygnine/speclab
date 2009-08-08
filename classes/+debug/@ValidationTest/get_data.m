function[data] = get_data(obj)
% DATA = VALIDATIONTEST.GET_DATA()
%
%     Runs the ValidationTest and sticks the boolean result in OBJ.RESULT.

data = obj.data_generator(obj.parameters);
