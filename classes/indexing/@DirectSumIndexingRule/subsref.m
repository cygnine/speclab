function[varargout] = subsref(self,s);
% subsref -- Subscripting for DirectSumIndexingRule
%     [varargout] = subsref(self,s);

if length(s)==1
  switch s(1).type
  case '()'
    [varargout{1:nargout}] = self.from_naturals(s.subs{:});
  case '.'
    try
      [varargout{1:nargout}] = self.(s(1).subs);
    catch 
      [varargout{1:nargout}] = builtin('subsref', self, s);
    end
  end
else
  % Call Singleton subsref
  [varargout{1:nargout}] = builtin('subsref', self, s);
end
