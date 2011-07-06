function[varargout] = subsref(self, s)
% subsref -- Overloaded subsref method for FactoredPolynomial class
%
% varargout = subsref(self, s)
%
%     s is a Matlab-standard subsref input. The purpose of this function is only
%     to output a help topic associated with a method if the method is called as
%     a script (i.e. without trailing parens). 

temp = methods(self);

% Then spit out help topic
if length(s)==1 
  if strcmpi(s(1).type, '.') && any(strcmpi(s(1).subs, temp))
    help(strcat(class(self), '/', s(1).subs));
  elseif strcmpi(s(1).type, '()')
    [varargout{1:nargout}] = self.evaluate(s.subs{:});
  else
    [varargout{1:nargout}] = builtin('subsref', self, s);
  end
else % Otherwise just do what you would normally do
  [varargout{1:nargout}] = builtin('subsref', self, s);
end
