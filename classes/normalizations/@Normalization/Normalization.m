classdef Normalization < Singleton
% Normalization -- A superclass for function normalizations 
%
% self = Normalization()
%
%     A superclass for all singleton classes that are normalization
%     specifications. There are no public data properties for this class, and
%     this class has no constructor.
%
% Normalization Properties:
%   ids - A set of strings and/or scalars that identify this normalization
% Normalization Methods:
%   string_compare - A method that tests if an input matches this normalization's ids
  properties(SetAccess=protected,Abstract=true)
    ids
  end
  properties(SetAccess=protected)
    description;
  end
  methods
    [bool, instance] = string_compare(self, string)
    bool = eq(self,other)
    disp(self)
  end
end
