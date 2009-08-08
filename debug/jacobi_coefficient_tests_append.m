function[container] = jacobi_coefficient_tests_append(container,opt);
% [CONTAINER] = JACOBI_COEFFICIENT_TESTS_APPEND(CONTAINER,OPT);
%
%     Appends some Jacobi coefficient ValidationTest's to the TestContainer
%     CONTAINER. The Jacobi family is defined by OPT. 

import debug.*

test = ValidationTest('description', '(1-r) x P coefficients',...
                      'parameters', opt,...
                      'validator', @one_minus_r_validator,...
                      'data_generator', @one_minus_r_data);
container = container.append(test);

test = ValidationTest('description', '(1+r) x P coefficients',...
                      'parameters', opt,...
                      'validator', @one_plus_r_validator,...
                      'data_generator', @one_plus_r_data);
container = container.append(test);

test = ValidationTest('description', '(1-r^2) x P coefficients',...
                      'parameters', opt,...
                      'validator', @one_minus_r_squared_validator,...
                      'data_generator', @one_minus_r_squared_data);
container = container.append(test);

if (opt.alpha>0) || (opt.beta>0)
  test = ValidationTest('description', 'Integer connection matrix',...
                        'parameters', opt,...
                        'validator', @integer_connection_validator,...
                        'data_generator', @integer_connection_data);
  container = container.append(test);
end

end

function[data] = one_minus_r_data(opt)
  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
  
  opt.N = min([opt.N,50]); % need not go above 50 polys
  promote_opt = opt;
  promote_opt.alpha = promote_opt.alpha+1;

  ns = 0:(promote_opt.N-1);
  ns_1 = 0:(promote_opt.N);
  jint = jac.interval(opt);
  r = linspace(jint(1),jint(2),300).';

  ps = jac.eval.eval_jacobi_poly(r,ns, promote_opt);
  mu = jac.coefficients.one_minus_r_times_p(ns,promote_opt.alpha,promote_opt.beta,promote_opt);
  ps_temp = jac.eval.eval_jacobi_poly(r,ns_1,opt);

  [data.ps,data.mu,data.ps_temp,data.r] = deal(ps,mu,ps_temp,r);
end

function[tf] = one_minus_r_validator(data,opt)
  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
  [ps,mu,ps_temp,r] = deal(data.ps, data.mu, data.ps_temp,data.r);
  opt.N = min([opt.N,50]); % need not go above 50 polys
  Nr = length(r);

  tol = 1e-5;
  tempopt = opt;
  tempopt.beta = 0;
  tempopt.alpha = 1;

  err = spdiags(jac.weights.weight(r,tempopt),0,Nr,Nr)*ps - ps_temp(:,1:opt.N)*spdiags(mu(:,1),0,opt.N,opt.N) - ...
                                ps_temp(:,2:(opt.N+1))*spdiags(mu(:,2),0,opt.N,opt.N);
  tf = max(max(abs(err)<tol));
end

function[data] = one_plus_r_data(opt)
  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
  
  opt.N = min([opt.N,50]); % need not go above 50 polys
  promote_opt = opt;
  promote_opt.beta = promote_opt.beta+1;

  ns = 0:(promote_opt.N-1);
  ns_1 = 0:(promote_opt.N);
  jint = jac.interval(opt);
  r = linspace(jint(1),jint(2),300).';

  ps = jac.eval.eval_jacobi_poly(r,ns, promote_opt);
  mu = jac.coefficients.one_plus_r_times_p(ns,promote_opt.alpha,promote_opt.beta,promote_opt);
  ps_temp = jac.eval.eval_jacobi_poly(r,ns_1,opt);

  [data.ps,data.mu,data.ps_temp,data.r] = deal(ps,mu,ps_temp,r);
end

function[tf] = one_plus_r_validator(data,opt)
  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
  [ps,mu,ps_temp,r] = deal(data.ps, data.mu, data.ps_temp,data.r);
  opt.N = min([opt.N,50]); % need not go above 50 polys
  Nr = length(r);

  tol = 1e-5;
  tempopt = opt;
  tempopt.beta = 1;
  tempopt.alpha = 0;

  err = spdiags(jac.weights.weight(r,tempopt),0,Nr,Nr)*ps - ps_temp(:,1:opt.N)*spdiags(mu(:,1),0,opt.N,opt.N) - ...
                                ps_temp(:,2:(opt.N+1))*spdiags(mu(:,2),0,opt.N,opt.N);

  tf = max(max(abs(err)<tol));
end

function[data] = one_minus_r_squared_data(opt)
  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
  
  opt.N = min([opt.N,50]); % need not go above 50 polys
  promote_opt = opt;
  promote_opt.alpha = promote_opt.alpha+1;
  promote_opt.beta = promote_opt.beta+1;

  ns = 0:(promote_opt.N-1);
  ns_2 = 0:(promote_opt.N+1);
  jint = jac.interval(opt);
  r = linspace(jint(1),jint(2),300).';

  ps = jac.eval.eval_jacobi_poly(r,ns, promote_opt);
  mu = jac.coefficients.one_minus_r_squared_times_p(ns,promote_opt.alpha,promote_opt.beta,promote_opt);
  ps_temp = jac.eval.eval_jacobi_poly(r,ns_2,opt);

  [data.ps,data.mu,data.ps_temp,data.r] = deal(ps,mu,ps_temp,r);
end

function[tf] = one_minus_r_squared_validator(data,opt)
  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
  [ps,mu,ps_temp,r] = deal(data.ps, data.mu, data.ps_temp,data.r);
  opt.N = min([opt.N,50]); % need not go above 50 polys
  Nr = length(r);

  tol = 1e-5;
  tempopt = opt;
  tempopt.beta = 1;
  tempopt.alpha = 1;

  err = spdiags(jac.weights.weight(r,tempopt),0,Nr,Nr)*ps - ps_temp(:,1:opt.N)*spdiags(mu(:,1),0,opt.N,opt.N) - ...
                                   ps_temp(:,2:(opt.N+1))*spdiags(mu(:,2),0,opt.N,opt.N) - ...
                                   ps_temp(:,3:end)*spdiags(mu(:,3),0,opt.N,opt.N); 

  tf = max(max(abs(err)<tol));
end

function[data] = integer_connection_data(opt)

  global handles;
  jac = handles.speclab.OrthogonalPolynomial1D.jacobi;

  f = @(x) sin(opt.N/(10*opt.scale)*x);

  [r,w] = jac.quad.gauss_quadrature(opt.N,opt);
  ps = jac.eval.eval_jacobi_poly(r,0:(opt.N-1), opt);
  gauss_modes = ps'*(w.*f(r));

  A = ceil(opt.alpha);
  B = ceil(opt.beta);
  demoted_opt = opt;
  demoted_opt.alpha = opt.alpha - A;
  demoted_opt.beta = opt.beta - B;
  [r_demoted,w_demoted] = jac.quad.gauss_quadrature(demoted_opt.N,demoted_opt);
  ps_demoted = jac.eval.eval_jacobi_poly(r_demoted,0:(opt.N-1),demoted_opt);
  modes = ps_demoted'*spdiags(w_demoted,0,opt.N,opt.N)*f(r_demoted);

  C = jac.connection.integer_separation_connection_matrix(opt.N,demoted_opt.alpha,demoted_opt.beta,A,B);
  [data.modes, data.C, data.demoted_opt, data.gauss_modes] = deal(modes,C,demoted_opt,gauss_modes);
end

function[tf] = integer_connection_validator(data,opt)
  
  [modes,C,demoted_opt,gauss_modes] = deal(data.modes, data.C, data.demoted_opt,data.gauss_modes);
  promoted_modes = C*modes;

  tol = 1e-5;
  AB = ceil(opt.alpha) + ceil(opt.beta);
  tf = max(abs(gauss_modes(1:(end-AB)) - promoted_modes(1:(end-AB))))<tol;
end
