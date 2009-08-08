function[] = extend(self,other)
% EXTEND(SELF, OTHER)
%
%     Add tests from TestContainer OTHER to tests of SELF.

if isa(other,'TestContainer')
  for n = 1:other.N
    self.append(other.tests{n});
  end
else
  fprintf('Error: input OTHER is not of type TestContainer\n');
end
