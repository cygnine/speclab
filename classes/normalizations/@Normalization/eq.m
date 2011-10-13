function[bool] = eq(self, other)
% eq -- Boolean equality comparison for FunctionNormalization 
%
% bool = eq(self, other)
%
%     Returns a true/false boolean by comparing self and other.
%
%     If other is a FunctionNormalization instance:
%        Returns true if the two instances are the same object.
%
%     If other is a string type:
%        Returns true if other matches any elements of self.ids

if isa(other, 'Normalization') && isa(self, 'Normalization');
  bool = eq@handle(self, other);
elseif isa(other, 'char')
  bool = any(strcmpi(other, self.ids));
elseif isa(self, 'char')
  bool = any(strcmpi(self, other.ids));
else
  error('Undefined method "==" for input type');
end
