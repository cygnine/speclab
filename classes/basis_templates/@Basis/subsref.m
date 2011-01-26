function[varargout] = subsref(self, s)
% subsref -- Overloaded subsref method for Basis class
%
% varargout = subsref(self, s)
%
%     s is a Matlab-standard subsref input. The purpose of this function is only
%     to output a help topic associated with a method if the method is called as
%     a script (i.e. without trailing parens). 

temp = methods(self);

% Then spit out help topic
if length(s)==1 && strcmp(s(1).type, '.') && any(strcmp(s(1).subs, temp))
  help(strcat(class(self), '/', s(1).subs));
else % Otherwise just do what you would normally do
  [varargout{1:nargout}] = builtin('subsref', self, s);
end
