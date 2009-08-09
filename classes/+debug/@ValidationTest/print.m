function[] = print(self,varargin)
% validationtest.print(self,{fmt=[]})
%
%     Pretty-prints the validation test description and parameters.

fprintf('%s\n', self.description);
fmt = varargin{1};

if ~strcmpi(fmt, 'brief')
  format compact
  fn = fieldnames(self.parameters);
  for n = 1:length(fn);
    fprintf('%15s : ', fn{n});
    disp(getfield(self.parameters, fn{n}));
  end
  format
end
