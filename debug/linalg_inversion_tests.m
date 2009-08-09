function[container] = linalg_inversion_tests(container,opt)
% [container] = linalg_inversion_tests(container,opt)
%
%     Runs la tests for inversion of operators.

import debug.*

test = ValidationTest('description', 'Inversion of sparse, upper triangular system',...
                      'parameters', opt,...
                      'validator', @sparse_triu_validator,...
                      'data_generator', @sparse_triu_data);
container = container.append(test);

end
function[data] = sparse_triu_data(opt)
  % If we allow negative values in s, it's horribly conditioned
  s = rand([opt.N,opt.bandwidth]);
  s = spdiags(s,0:(opt.bandwidth-1),opt.N,opt.N);
  b = randn([opt.N,1]);
  [data.s, data.b] = deal(s,b);
end
  
function[tf] = sparse_triu_validator(data,opt)
  global handles;
  linv = handles.common.linalg.triu_sparse_invert;

  x1 = linv(data.s, data.b, 'bandwidth', opt.bandwidth);
  x2 = inv(data.s)*data.b;

  tol = 10^(-10 + opt.bandwidth);

  tf = max(abs(x1-x2))<tol;
end
