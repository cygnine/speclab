function self = extend(self,other)
% EXTEND(SELF, OTHER)
%
%     Add tests from TestContainer OTHER to tests of SELF.

import debug.TestContainer

if isa(other,'TestContainer')
  for n = 1:other.N
    self = self.append(other.tests{n});
  end
else
  fprintf('Error: input OTHER is not of type TestContainer\n');
end
