function[container] = laguerre_coefficient_tests_append(container,opt);
% [container] = laguerre_coefficient_tests_append(container,opt);
%
%     Appends some laguerre coefficient ValidationTest's to the TestContainer
%     container. The laguerre family is defined by opt. 

import debug.*

test = ValidationTest('description', 'x x p coefficients',...
                      'parameters', opt,...
                      'validator', @x_times_p_validator,...
                      'data_generator', @x_times_p_data);
container = container.append(test);

if (opt.alpha>0)
  test = ValidationTest('description', 'Integer connection matrix',...
                        'parameters', opt,...
                        'validator', @integer_connection_validator,...
                        'data_generator', @integer_connection_data);
  container = container.append(test);
end

test = ValidationTest('description', 'd/dr P coefficients',...
                      'parameters', opt,...
                      'validator', @ddr_P_validator,...
                      'data_generator', @ddr_P_data);
container = container.append(test);

function[data] = x_times_p_data(opt)
  global handles;
  lag = handles.speclab.orthopoly1d.laguerre;
  
  opt.N = min([opt.N,50]); % need not go above 50 polys
  promote_opt = opt;
  promote_opt.alpha = promote_opt.alpha+1;

  ns = 0:(promote_opt.N-1);
  ns_1 = 0:(promote_opt.N);
  lint = lag.interval(opt);
  if lint(2)==Inf
    r = linspace(opt.shift,opt.shift+2*opt.scale,300).';
  else
    r = linspace(opt.shift-2*opt.scale, opt.shift, 300).';
  end

  ps = lag.eval.eval_laguerre_poly(r,ns, promote_opt);
  mu = lag.coefficients.x_times_p(ns,promote_opt.alpha,promote_opt);
  ps_temp = lag.eval.eval_laguerre_poly(r,ns_1,opt);

  [data.ps,data.mu,data.ps_temp,data.r] = deal(ps,mu,ps_temp,r);

function[tf] = x_times_p_validator(data,opt)
  global handles;
  lag = handles.speclab.orthopoly1d.laguerre;
  sss = handles.speclab.common.standard_scaleshift_1d;
  [ps,mu,ps_temp,r] = deal(data.ps, data.mu, data.ps_temp,data.r);
  opt.N = min([opt.N,50]); % need not go above 50 polys
  Nr = length(r);

  tol = 1e-5;
  tempopt = opt;
  tempopt.alpha = 1;

  r = sss(r,opt);
  err = spdiags(r,0,Nr,Nr)*ps - ps_temp(:,1:opt.N)*spdiags(mu(:,1),0,opt.N,opt.N) - ...
                                ps_temp(:,2:(opt.N+1))*spdiags(mu(:,2),0,opt.N,opt.N);
  tf = max(max(abs(err)<tol));

function[data] = integer_connection_data(opt)

  global handles;
  lag = handles.speclab.orthopoly1d.laguerre;

  f = @(x) exp(-x.^2);

  [r,w] = lag.quad.gauss_quadrature(opt.N,opt);
  ps = lag.eval.eval_laguerre_poly(r,0:(opt.N-1), opt);
  gauss_modes = ps'*(w.*f(r));

  A = ceil(opt.alpha);
  demoted_opt = opt;
  demoted_opt.alpha = opt.alpha - A;
  [r_demoted,w_demoted] = lag.quad.gauss_quadrature(demoted_opt.N,demoted_opt);
  ps_demoted = lag.eval.eval_laguerre_poly(r_demoted,0:(opt.N-1),demoted_opt);
  modes = ps_demoted'*spdiags(w_demoted,0,opt.N,opt.N)*f(r_demoted);

  C = lag.connection.integer_separation_connection_matrix(opt.N,demoted_opt.alpha,A);
  [data.modes, data.C, data.demoted_opt, data.gauss_modes] = deal(modes,C,demoted_opt,gauss_modes);

function[tf] = integer_connection_validator(data,opt)
  
  [modes,C,demoted_opt,gauss_modes] = deal(data.modes, data.C, data.demoted_opt,data.gauss_modes);
  promoted_modes = C*modes;

  tol = 1e-3;
  AB = ceil(opt.alpha);
  tf = max(abs(gauss_modes(1:(end-AB)) - promoted_modes(1:(end-AB))))<tol;

function[data] = ddr_P_data(opt);
  global handles;
  lag = handles.speclab.orthopoly1d.laguerre;

  lint = lag.interval(opt);
  if lint(2)==Inf
    r = linspace(opt.shift, opt.shift+2*opt.scale,300).';
  else
    r = linspace(opt.shift-2*opt.scale, opt.shift, 300).'
  end
  zetas = lag.coefficients.derivative(0:(opt.N-1), opt.alpha, opt);
  opt.d = 1;
  dps = lag.eval.eval_laguerre_poly(r,0:(opt.N-1),opt);
  opt.d = 0;
  opt.alpha = opt.alpha+1;
  ps = lag.eval.eval_laguerre_poly(r,0:(opt.N-1),opt);
  [data.r, data.ps, data.dps, data.zetas] = deal(r,ps,dps,zetas);

function[tf] = ddr_P_validator(data,opt);
  [r,ps,dps,zetas] = deal(data.r, data.ps, data.dps, data.zetas);

  tol = 10^(-8+(opt.alpha/3));
  opt.N = min([opt.N, 50]); % don't test them all
  err = dps(:,2:opt.N) - ps(:,1:(opt.N-1))*spdiags(zetas(2:end),0,opt.N-1,opt.N-1);
  tf = norm(err)<tol;
