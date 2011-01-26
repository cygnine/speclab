function[varargout] = subsref(self,s);
% subsref -- Subsref method for AffineMap
%
% varargout = subsref(self,s)

if length(s)==1
  switch s(1).type
  case '()'
    if isa(s.subs{1}, 'AffineMap');
      [varargout{1:nargout}] = self.compose(s.subs{:});
    else
      [varargout{1:nargout}] = self.evaluate(s.subs{:});
    end
  case '.'
    switch s(1).subs
    %case {'m', 'n'}
    %  error('Cannot access private properties')
    otherwise
      [varargout{1:nargout}] = self.(s(1).subs);
    end
  otherwise 
    error('Unrecognized subscripting of function');
  end
%elseif length(s)==2  % This must be self.handle(input1, input2, ...)
%  [varargout{1:nargout}] = self.handle(s(2).subs{:});
elseif length(s)==2
  if strcmp(s(1).type, '.') & strcmp(s(2).type, '()')
    [varargout{1:nargout}] = self.(s(1).subs)(s(2).subs{:});
  else
    fprintf('Not implemented for multiple-subscripting of AffineMap\n');
  end
else
  fprintf('Not implemented for multiple-subscripting of AffineMap\n');
end
