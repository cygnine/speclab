function[container] = jacobi_operator_tests_append(container,opt);
% [container] = jacobi_operator_tests_append(container,opt);
%
%     Appends some Jacobi coefficient ValidationTest's to the TestContainer
%     container. The Jacobi family is defined by opt. 

import debug.*

test = ValidationTest('description', 'Stiffness operator',...
                      'parameters', opt,...
                      'validator', @stiffness_validator,...
                      'data_generator', @stiffness_data);
container = container.append(test);

end

function[data] = stiffness_data(opt);
  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;

  modes = randn([opt.N,1]);
  [r,w] = jac.quad.gauss_quadrature(opt.N,opt);
  opt.d = 1;
  dps = jac.eval.eval_jacobi_poly(r,0:(opt.N-1),opt);
  opt.d = 0;
  ps = jac.eval.eval_jacobi_poly(r,0:(opt.N-1),opt);
  vinv = ps'*spdiags(w,0,opt.N,opt.N);

  dfun_rec = dps*modes;
  dmodes = vinv*dfun_rec;

  [data.modes,data.dmodes] = deal(modes,dmodes);
end

function[tf] = stiffness_validator(data,opt);
  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;

  [modes,dmodes] = deal(data.modes,data.dmodes);

  dmodes2 = jac.operators.stiffness_operator(modes,opt);
  tol = 10^(-8 + opt.alpha/3 + opt.beta/3 + abs(opt.alpha-opt.beta)/3);
  tf = norm(dmodes - dmodes2)<tol;
end
