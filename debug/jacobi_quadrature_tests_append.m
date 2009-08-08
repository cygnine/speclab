function[container] = jacobi_quadrature_tests_append(container,opt)
% [CONTAINER] = JACOBI_QUADRATURE_TESTS_APPEND(CONTAINER,OPT)
%
%     Called by jacobi_tests to add ValidationTests given a base Jacobi
%     polynomial family descriptor OPT.

import debug.*

test = ValidationTest('description', 'Gauss quadrature nodes inside interval',...
                      'parameters', opt,...
                      'validator', @gauss_quad_nodes,...
                      'data_generator', @gauss_quad_data);
container = container.append(test);

test = ValidationTest('description', 'Gauss quadrature weights positive',...
                      'parameters', opt,...
                      'validator', @gauss_quad_weights,...
                      'data_generator', @gauss_quad_data);
container = container.append(test);

opt.r = -opt.scale+opt.shift;
opt.r1 = opt.r;
test = ValidationTest('description', 'Gauss-Radau quadrature nodes inside interval',...
                      'parameters', opt,...
                      'validator', @gauss_quad_nodes,...
                      'data_generator', @gauss_radau_data);
container = container.append(test);

test = ValidationTest('description', 'Gauss-Radau quadrature weights positive',...
                      'parameters', opt,...
                      'validator', @gauss_quad_weights,...
                      'data_generator', @gauss_radau_data);
container = container.append(test);

opt.r = opt.scale+opt.shift;
opt.r2 = opt.r;
test = ValidationTest('description', 'Gauss-Radau quadrature nodes inside interval',...
                      'parameters', opt,...
                      'validator', @gauss_quad_nodes,...
                      'data_generator', @gauss_radau_data);
container = container.append(test);

test = ValidationTest('description', 'Gauss-Radau quadrature weights positive',...
                      'parameters', opt,...
                      'validator', @gauss_quad_weights,...
                      'data_generator', @gauss_radau_data);
container = container.append(test);

test = ValidationTest('description', 'Gauss-Lobatto quadrature nodes inside interval',...
                      'parameters', opt,...
                      'validator', @gauss_quad_nodes,...
                      'data_generator', @gauss_radau_data);
container = container.append(test);

test = ValidationTest('description', 'Gauss-Lobatto quadrature weights positive',...
                      'parameters', opt,...
                      'validator', @gauss_quad_weights,...
                      'data_generator', @gauss_lobatto_data);
container = container.append(test);

end

function[data] = gauss_quad_data(opt)
  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;

  [x,w] = jac.quad.gauss_quadrature(opt.N,opt);
  data.x = x;
  data.w = w;
end

function[data] = gauss_radau_data(opt)
  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;

  [x,w] = jac.quad.gauss_radau_quadrature(opt.N,opt);
  data.x = x;
  data.w = w;
end

function[data] = gauss_lobatto_data(opt)
  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;

  [x,w] = jac.quad.gauss_lobatto_quadrature(opt.N,opt);
  data.x = x;
  data.w = w;
end

function[tf] = gauss_quad_nodes(data,opt)
  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;

  jint = jac.interval(opt);
  tol = 1e-12;
  [x,w] = deal(data.x,data.w);
  tf = all(x<jint(2)+tol) && all(x>jint(1)-tol);
end

function[tf] = gauss_quad_weights(data,opt)
  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;

  tol = 1e-12;
  jint = jac.interval(opt);
  [x,w] = deal(data.x,data.w);
  tf = all(w>-tol);
end
