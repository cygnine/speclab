function[container] = wiener_matrix_tests(container,opt)
% [container] = wiener_matrix_tests(container,opt)
%
%     Appends function matrix tests for Wiener functions to the TestContainer
%     container. The Wiener family to use is specified by the struct of
%     parameters opt.

import debug.*

test = ValidationTest('description', 'Wiener stiffness matrix',...
                      'parameters', opt,...
                      'validator', @stiff_validator,...
                      'data_generator', @stiff_data);
container = container.append(test);

function[data] = stiff_data(opt)
  global packages;
  wiener = packages.speclab.wiener;

  [x,w] = wiener.quad.pi_gauss_quadrature(2*opt.N,opt);
  ks = packages.speclab.common.integer_range(opt.N);
  ws = wiener.eval.wiener_function(x,ks,opt);
  dws = wiener.eval.derivative_wiener_function(x,ks,opt);

  stiff_quadrature = ws'*spdiags(w,0,2*opt.N,2*opt.N)*dws;
  stiff_sparse = wiener.matrices.wiener_stiffness_matrix(opt.N,opt);
  [data.stiff_quadrature, data.stiff_sparse] = ...
     deal(stiff_quadrature, stiff_sparse);

function[tf] = stiff_validator(data,opt)
  [stiff_quadrature, stiff_sparse] = ...
    deal(data.stiff_quadrature, data.stiff_sparse);

  tol = 1e-6;
  tf = norm(stiff_quadrature - stiff_sparse)<tol;
