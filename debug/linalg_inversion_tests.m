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

function[data] = sparse_triu_data(opt)
  % If we allow negative values in s, it's horribly conditioned
  s = rand([opt.N,opt.bandwidth]);
  s = spdiags(s,0:(opt.bandwidth-1),opt.N,opt.N);
  b = randn([opt.N,1]);
  [data.s, data.b] = deal(s,b);
  
function[tf] = sparse_triu_validator(data,opt)
  global packages;
  linv = packages.labtools.linalg.triu_sparse_invert;

  x1 = linv(data.s, data.b, 'bandwidth', opt.bandwidth);
  x2 = inv(data.s)*data.b;

  tol = min([10^(-8 + opt.bandwidth), 0.1]);

  tf = max(abs(x1-x2))<tol;
