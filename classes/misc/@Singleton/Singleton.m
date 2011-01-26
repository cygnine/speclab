classdef Singleton < handle;
% Singleton -- Defines a Singleton class (a class wtih only one instance)
%
% self = Singleton()
%
%     A boilerplate class definition for a singleton object -- that is, an
%     object that is instantiated once and only once. This is meant as a
%     superclass for other custom singleton classes.
%
%     This implementation is based upon the "Singleton" class implementation
%     provided by The Mathworks. This original implementation can be accessed at 
%     http://www.mathworks.com.au/matlabcentral/fileexchange/24911.

  methods(Abstract, Static)
    self = instance();
  end
end
