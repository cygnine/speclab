function[container] = jacobi_approximation_tests_append(container,opt);
% [CONTAINER] = JACOBI_APPROXIMATION_TESTS_APPEND(CONTAINER,OPT);
%
%     Appends some Jacobi approximation ValidationTest's to the TestContainer
%     CONTAINER. The Jacobi family is defined by OPT. 

import debug.*

opt.n = 0:(opt.N-1);
test = ValidationTest('description', 'Mass matrix identity, Gauss quadrature',...
                      'parameters', opt,...
                      'validator', @mass_matrix_validator,...
                      'data_generator', @mass_data);
container = container.append(test);

test = ValidationTest('description', 'Jacobi Interpolant approximation',...
                      'parameters', opt,...
                      'validator', @interpolant_validator,...
                      'data_generator', @interpolant_data);
container = container.append(test);

test = ValidationTest('description', 'Jacobi Derivative approximation',...
                      'parameters', opt,...
                      'validator', @derivative_validator,...
                      'data_generator', @derivative_data);
container = container.append(test);

test = ValidationTest('description', 'Jacobi Derivative approximation (stiffness)',...
                      'parameters', opt,...
                      'validator', @derivative_stiffness_validator,...
                      'data_generator', @derivative_stiffness_data);
container = container.append(test);

function[data] = mass_data(opt)
  global packages;
  jac = packages.speclab.orthopoly1d.jacobi;

  [x,w] = jac.quad.gauss_quadrature(opt.N,opt);
  ps = jac.eval.eval_jacobi_poly(x,opt.n,opt);
  [data.x,data.w,data.ps] = deal(x,w,ps);

function[data] = interpolant_data(opt)
  global packages;
  jac = packages.speclab.orthopoly1d.jacobi;
  jint = jac.interval(opt);

  [x,w] = jac.quad.gauss_quadrature(opt.N,opt);
  ps = jac.eval.eval_jacobi_poly(x,opt.n,opt);

  x_refined = linspace(jint(1),jint(2),5*opt.N)';
  ps_refined = jac.eval.eval_jacobi_poly(x_refined,opt.n,opt);
  [data.x,data.w,data.ps,data.x_refined,data.ps_refined] = deal(...
    x,w,ps,x_refined,ps_refined);

function[data] = derivative_data(opt)
  global packages;
  jac = packages.speclab.orthopoly1d.jacobi;
  jint = jac.interval(opt);

  [x,w] = jac.quad.gauss_quadrature(opt.N,opt);
  ps = jac.eval.eval_jacobi_poly(x,opt.n,opt);

  x_refined = linspace(jint(1),jint(2),5*opt.N)';
  x_refined = x_refined(2:end-1);
  opt.d = 1;
  dps_refined = jac.eval.eval_jacobi_poly(x_refined,opt.n,opt);
  [data.x,data.w,data.ps,data.x_refined,data.dps_refined] = deal(...
    x,w,ps,x_refined,dps_refined);

function[tf] = mass_matrix_validator(data,opt)
  [x,w,ps] = deal(data.x,data.w,data.ps);

  tol = 10^(-8+abs(opt.alpha-opt.beta)/5 + (opt.alpha+opt.beta)/8);
  mass = ps'*spdiags(w,0,opt.N,opt.N)*ps;
  tf = norm(mass-eye(opt.N))<tol;

function[tf] = interpolant_validator(data,opt);

  f = @(x) sin(opt.N/(10*opt.scale)*x);
  tol = 10^(-7+opt.alpha/3+opt.beta/3 + abs(opt.alpha-opt.beta)/3);
  [x,w,ps,x_refined,ps_refined] = deal(data.x,data.w,data.ps, ...
     data.x_refined, data.ps_refined);

  modes = ps'*(w.*f(x));
  f_approx = ps_refined*modes;
  tf = all(abs(f_approx - f(x_refined))<tol);

function[tf] = derivative_validator(data,opt);

  global packages
  jac = packages.speclab.orthopoly1d.jacobi;
  f = @(x) sin(opt.N/(10*opt.scale)*x);
  df = @(x) opt.N/(10*opt.scale)*cos(opt.N/(10*opt.scale)*x);
  tol = 2*10^(-4+opt.alpha/6+opt.beta/6 + abs(opt.alpha-opt.beta)/6);
  [x,w,ps,x_refined,dps_refined] = deal(data.x,data.w,data.ps, ...
     data.x_refined, data.dps_refined);

  modes = ps'*(w.*f(x));
  df_approx = dps_refined*modes;
  wgt = sqrt(jac.weights.weight(x_refined,opt));
  tf = all(abs(df_approx - df(x_refined)).*wgt<tol);

function[data] = derivative_stiffness_data(opt);
  global packages;
  jac = packages.speclab.orthopoly1d.jacobi;
  jint = jac.interval(opt);

  [x,w] = jac.quad.gauss_quadrature(opt.N,opt);
  ps = jac.eval.eval_jacobi_poly(x,opt.n,opt);

  vinv = ps'*spdiags(w,0,opt.N,opt.N);

  x_refined = linspace(jint(1),jint(2),5*opt.N)';
  [data.x,data.x_refined,data.vinv,data.ps] = deal(...
    x,x_refined,vinv,ps);

function[tf] = derivative_stiffness_validator(data,opt);
  global packages;
  jac = packages.speclab.orthopoly1d.jacobi;

  f = @(x) sin(opt.N/(10*opt.scale)*x);
  df = @(x) opt.N/(10*opt.scale)*cos(opt.N/(10*opt.scale)*x);
  tol = 9*10^(-5+opt.alpha/6+opt.beta/6 + abs(opt.alpha-opt.beta)/6);
  [x,x_refined,vinv,ps] = deal(data.x,data.x_refined,data.vinv,data.ps);
  modes = vinv*f(x);

  dmodes = jac.operators.stiffness_operator(modes,opt);
  drec = ps*dmodes;

  tf = norm(drec-df(x))<tol;
