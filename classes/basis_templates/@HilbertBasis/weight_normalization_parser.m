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
bool = false;
q = 1;
while not(bool) & (q <= length(self.allowed_weight_normalizations))
  curr_norm = self.allowed_weight_normalizations{q};
  [bool, obj] = curr_norm.string_compare(inp);
  q = q + 1;
end
if not(bool)
  obj = self.default_weight_normalization;
  %error('Unrecognized normalization character array');
end

%switch lower(inp)
%case {'classical', 'class'}
%  obj = ClassicalWeightNormalization.instance();
%end
