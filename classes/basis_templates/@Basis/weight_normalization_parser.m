function[obj] = weight_normalization_parser(self, inp)
% weight_normalization_parser -- Returns a WeightNormalization subclass
%
% obj = weight_normalization_parser(input)
% 
%     The given input is a string (or a WeightNormalization subclass) that
%     specifies some kind of known weight normalization type. This function
%     parses the input and returns the appropriate WeightNormalization
%     subclass instance, if appropriate. 
%
%     If the input type is not recognized, an error is NOT thrown -- instead, an
%     empty array [] is returned.

if not(isa(inp, 'char'))
  if isa(inp, 'WeightNormalization')
    obj = inp;
  else
    obj = [];
  end
  return
end

obj = [];
switch lower(inp)
case {'classical', 'class'}
  obj = ClassicalWeightNormalization.instance();
end