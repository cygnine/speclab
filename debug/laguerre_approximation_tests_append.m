function[container] = laguerre_approximation_tests_append(container,opt);
% [container] = laguerre_approximation_tests_append(container,opt);
%
%     Appends some laguerre approximation ValidationTest's to the TestContainer
%     container. The laguerre family is defined by opt. 

import debug.*

opt.n = 0:(opt.N-1);
test = ValidationTest('description', 'Mass matrix identity, Gauss quadrature',...
                      'parameters', opt,...
                      'validator', @mass_matrix_validator,...
                      'data_generator', @mass_data);
container = container.append(test);

test = ValidationTest('description', 'laguerre Interpolant approximation',...
                      'parameters', opt,...
                      'validator', @interpolant_validator,...
                      'data_generator', @interpolant_data);
container = container.append(test);

test = ValidationTest('description', 'laguerre Derivative approximation',...
                      'parameters', opt,...
                      'validator', @derivative_validator,...
                      'data_generator', @derivative_data);
container = container.append(test);

test = ValidationTest('description', 'laguerre Derivative approximation (stiffness)',...
                      'parameters', opt,...
                      'validator', @derivative_stiffness_validator,...
                      'data_generator', @derivative_stiffness_data);
container = container.append(test);

function[data] = mass_data(opt)
  from speclab.orthopoly1d import laguerre as lag
  %lag = packages.speclab.orthopoly1d.laguerre;

  [x,w] = lag.quad.gauss_quadrature(opt.N,opt);
  ps = lag.eval.eval_laguerre_poly(x,opt.n,opt);
  [data.x,data.w,data.ps] = deal(x,w,ps);

function[data] = interpolant_data(opt)
  from speclab.orthopoly1d import lageurre as lag
  %lag = packages.speclab.orthopoly1d.laguerre;
  lint = lag.interval(opt);

  [x,w] = lag.quad.gauss_quadrature(opt.N,opt);
  ps = lag.eval.eval_laguerre_poly(x,opt.n,opt);

  if lint(2)==Inf
    %x_refined = linspace(lint(1), x(end), 5*opt.N).';
    x_refined = linspace(lint(1), opt.scale*10, 5*opt.N).';
  else
    x_refined = linspace(x(1), lint(2), 5*opt.N).';
  end
  ps_refined = lag.eval.eval_laguerre_poly(x_refined,opt.n,opt);
  [data.x,data.w,data.ps,data.x_refined,data.ps_refined] = deal(...
    x,w,ps,x_refined,ps_refined);

function[data] = derivative_data(opt)
  from speclab.orthopoly1d import lageurre as lag
  %lag = packages.speclab.orthopoly1d.laguerre;
  lint = lag.interval(opt);

  [x,w] = lag.quad.gauss_quadrature(opt.N,opt);
  ps = lag.eval.eval_laguerre_poly(x,opt.n,opt);

  if lint(2)==Inf
    %x_refined = linspace(lint(1), x(end), 5*opt.N).';
    x_refined = linspace(lint(1), opt.scale*10, 5*opt.N).';
  else
    x_refined = linspace(x(1), lint(2), 5*opt.N).';
  end
  opt.d = 1;
  dps_refined = lag.eval.eval_laguerre_poly(x_refined,opt.n,opt);
  [data.x,data.w,data.ps,data.x_refined,data.dps_refined] = deal(...
    x,w,ps,x_refined,dps_refined);

function[tf] = mass_matrix_validator(data,opt)
  [x,w,ps] = deal(data.x,data.w,data.ps);

  tol = 10^(-8+abs(opt.alpha)/5 + (opt.alpha)/8);
  mass = ps'*spdiags(w,0,opt.N,opt.N)*ps;
  tf = norm(mass-eye(opt.N))<tol;

function[tf] = interpolant_validator(data,opt);
  from speclab.orthopoly1d import lageurre as lag
  %lag = packages.speclab.orthopoly1d.laguerre;

  f = @(x) exp(-x.^2);
  tol = 10^(-3.5+opt.alpha/3);
  [x,w,ps,x_refined,ps_refined] = deal(data.x,data.w,data.ps, ...
     data.x_refined, data.ps_refined);

  modes = ps'*(w.*f(x));
  f_approx = ps_refined*modes;
  wgt = sqrt(lag.weights.weight(x_refined,opt));
  tf = all(abs(f_approx - f(x_refined)).*wgt<tol);

function[tf] = derivative_validator(data,opt);

  from speclab.orthopoly1d import lageurre as lag
  %lag = packages.speclab.orthopoly1d.laguerre;
  f = @(x) exp(-x.^2);
  df = @(x) -2*x.*exp(-x.^2);
  tol = 2*10^(-2+opt.alpha/6);
  [x,w,ps,x_refined,dps_refined] = deal(data.x,data.w,data.ps, ...
     data.x_refined, data.dps_refined);

  modes = ps'*(w.*f(x));
  df_approx = dps_refined*modes;
  wgt = sqrt(lag.weights.weight(x_refined,opt));
  tf = all(abs(df_approx - df(x_refined)).*wgt<tol);

function[data] = derivative_stiffness_data(opt);
  from speclab.orthopoly1d import lageurre as lag
  %lag = packages.speclab.orthopoly1d.laguerre;
  jint = lag.interval(opt);

  [x,w] = lag.quad.gauss_quadrature(opt.N,opt);
  ps = lag.eval.eval_laguerre_poly(x,opt.n,opt);

  vinv = ps'*spdiags(w,0,opt.N,opt.N);

  x_refined = linspace(jint(1),jint(2),5*opt.N)';
  [data.x,data.x_refined,data.vinv,data.ps] = deal(...
    x,x_refined,vinv,ps);

function[tf] = derivative_stiffness_validator(data,opt);
  from speclab.orthopoly1d import lageurre as lag
  %lag = packages.speclab.orthopoly1d.laguerre;

  f = @(x) exp(-x.^2);
  df = @(x) -2*x.*exp(-x.^2);
  tol = 2*10^(-2+opt.alpha/6);
  [x,x_refined,vinv,ps] = deal(data.x,data.x_refined,data.vinv,data.ps);
  modes = vinv*f(x);

  dmodes = lag.operators.stiffness_operator(modes,opt);
  drec = ps*dmodes;

  wgt = sqrt(lag.weights.weight(x,opt));
  tf = norm((drec-df(x)).*wgt)<tol;
