classdef Normalization < Singleton
% Normalization -- A superclass for function normalizations 
%
% self = Normalization()
%
%     A superclass for all singleton classes that are normalization
%     specifications. There are no public data properties for this class, and
%     this class has no constructor.
  properties(SetAccess=private,Abstract=true)
    ids
  end
  methods
    [bool, instance] = string_compare(self, string)
  end
end
