% Tests the connection coefficients mu

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
Nq = 500;
N = 300;
ns = (0:(N-1)).';
ns_plus_1 = (0:N).';
ns_plus_2 = (0:(N+1)).';

f = @(x) sin(100*x);
df = @(x) 100*cos(100*x);

alpha = 30*rand;
alpha = -1/2;
beta = 30*rand;
beta = 3/2;
opt.alpha = alpha;
opt.beta = beta;

[r,w] = jac.gauss_quadrature(Nq,opt);
ps = jac.eval_jacobi_poly(r,ns,opt);

% Test (1-r) x p:
if opt.alpha > 0
  mu = jac.coefficients.one_minus_r_times_p(ns,opt.alpha,opt.beta);

  ps_temp = jac.eval_jacobi_poly(r,ns_plus_1,'alpha',opt.alpha-1, ...
                                             'beta', opt.beta);

  err = spdiags(1-r,0,Nq,Nq)*ps - ps_temp(:,1:N)*spdiags(mu(:,1),0,N,N) - ...
                                  ps_temp(:,2:(N+1))*spdiags(mu(:,2),0,N,N);

  fprintf('Error for (1-r) x p is %1.4e\n', norm(err));
end

% Test (1+r) x p:
if opt.beta >0
  mu = jac.coefficients.one_plus_r_times_p(ns,opt.alpha,opt.beta);
  ps_temp = jac.eval_jacobi_poly(r,ns_plus_1,'alpha', opt.alpha,...
                                             'beta', opt.beta-1);

  err = spdiags(1+r,0,Nq,Nq)*ps - ps_temp(:,1:N)*spdiags(mu(:,1),0,N,N) - ...
                                  ps_temp(:,2:(N+1))*spdiags(mu(:,2),0,N,N);

  fprintf('Error for (1+r) x p is %1.4e\n', norm(err));
end

% Test (1-r^2) x p:
if (opt.alpha>0) & (opt.beta>0)
  mu = jac.coefficients.one_minus_r_squared_times_p(ns,opt.alpha,opt.beta);
  ps_temp = jac.eval_jacobi_poly(r,ns_plus_2,'alpha', opt.alpha-1,...
                                             'beta', opt.beta-1);

  err = spdiags(1-r.^2,0,Nq,Nq)*ps - ps_temp(:,1:N)*spdiags(mu(:,1),0,N,N) - ...
                                     ps_temp(:,2:(N+1))*spdiags(mu(:,2),0,N,N) - ...
                                     ps_temp(:,3:end)*spdiags(mu(:,3),0,N,N); 

  fprintf('Error for (1-r^2) x p is %1.4e\n', norm(err));
end

% Test modal connections:
gauss_modes = ps'*spdiags(w,0,Nq,Nq)*f(r);

A = ceil(opt.alpha);
B = ceil(opt.beta);
alpha_demoted = opt.alpha - A;
beta_demoted = opt.beta - B;
[r_demoted,w_demoted] = jac.gauss_quadrature(Nq,'alpha',alpha_demoted,'beta',beta_demoted);
ps_demoted = jac.eval_jacobi_poly(r_demoted,ns,'alpha',alpha_demoted,'beta',beta_demoted);
modes = ps_demoted'*spdiags(w_demoted,0,Nq,Nq)*f(r_demoted);

C = jac.connection.integer_separation_connection_matrix(N,alpha_demoted,beta_demoted,A,B);
promoted_modes = C*modes;

fprintf('Error for modal promotion is %1.4e\n', norm(gauss_modes -...
  promoted_modes));
