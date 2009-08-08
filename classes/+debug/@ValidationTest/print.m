function[] = print(obj)
% VALIDATIONTEST.PRINT()
%
%     Pretty-prints the validation test description and parameters.

fprintf('%s\n', obj.description);

fn = fieldnames(obj.parameters);
format compact
for n = 1:length(fn);
  fprintf('    %15s : ', fn{n});
  disp(getfield(obj.parameters, fn{n}));
end
format
