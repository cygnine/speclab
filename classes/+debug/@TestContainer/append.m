function self = append(self,test)
% SELF = APPEND(SELF,TEST)
%
%     Adds input TEST (type ValidationTest) to the container SELF.

import debug.ValidationTest

if isa(test,'ValidationTest')
  self.tests{end+1} = test;
  self.N = self.N + 1;
else
  fprintf('Error: input TEST is not of type ValidationTest\n')
end
