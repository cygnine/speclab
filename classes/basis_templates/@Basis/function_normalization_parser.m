function[obj] = function_normalization_parser(self,inp)
% function_normalization_parser -- Returns a FunctionNormalization subclass
%
% obj = function_normalization_parser(input)
% 
%     The given input is a string (or a FunctionNormalization subclass) that
%     specifies some kind of known function normalization type. This function
%     parses the input and returns the appropriate FunctionNormalization
%     subclass instance, if appropriate. 
%
%     If the input type is not recognized, an error is NOT thrown -- instead, an
%     empty array [] is returned.

if not(isa(inp, 'char'))
  if isa(inp, 'FunctionNormalization')
    obj = inp;
  else
    error('I don''t know what to do with this non-character function normalization');
  end
  return
end

obj = [];
bool = false;
q = 1;
while not(bool) & (q <= length(self.allowed_function_normalizations))
  curr_norm = self.allowed_function_normalizations{q};
  [bool, obj] = curr_norm.string_compare(inp);
  q = q + 1;
end
if not(bool)
  obj = self.default_function_normalization;
  warning(['I don''t know the normalization ID "' inp '", but I''m ignoring it']);
  %error('Unrecognized normalization character array');
end

%switch lower(inp)
%%case {'classical', 'class'}
%%  obj = ClassicalFunctionNormalization.instance();
%%case {'orthonormal', 'normal'}
%%  obj = OrthonormalNormalization.instance();
%%case {'monic', 'mon'}
%%  obj = MonicNormalization.instance();
%otherwise
%end
