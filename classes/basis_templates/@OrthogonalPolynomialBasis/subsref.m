function[varargout] = subsref(self,s);
% subsref -- Subscripting for OrthogonalPolynomialBasis
%     [varargout] = subsref(self,s);

if length(s)==1
  switch s(1).type
  case '()'
    [varargout{1:nargout}] = self.evaluate(s.subs{:});
  case '.'
    try
      [varargout{1:nargout}] = self.(s(1).subs);
    catch 
      [varargout{1:nargout}] = subsref@Basis(self, s);
    end
  otherwise 
    %error('Unrecognized subscripting of function');
    [varargout{1:nargout}] = subsref@Basis(self, s);
  end
elseif (length(s)==2) && strcmp(s(1).type, '.') & strcmp(s(2).type, '()')
  [varargout{1:nargout}] = self.(s(1).subs)(s(2).subs{:});
else
  % Call Basis subsref
  [varargout{1:nargout}] = subsref@Basis(self, s);
end
