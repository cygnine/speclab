function[bool, instance] = string_compare(self, string)
% string_compare -- Compares input string to instance 'ids' property
%
% [bool, instance] = string_compare(self, string)
%
%     Given a string 'string', this function sets bool=true if 'string' is a
%     case-insensitive match to any of the cell entries of self.ids.

bool = self==string;
if bool
  instance = self;
else
  instance = [];
end
