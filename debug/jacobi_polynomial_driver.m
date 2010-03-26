function[varargout] = jacobi_polynomial_driver(varargin)
% [FLAGS,DESCRIPTIONS] = JACOBI_POLYNOMIAL_DRIVER({N=10,ALPHA=-1/2,BETA=-1/2,SCALE=1,SHIFT=0})
% 
%     Runs various tests on speclab's Jacobi polynomial package. FLAGS is a
%     boolean array containing success indicators for each test, whose
%     description is in the cell array DESCRIPTIONS.

global packages;
jac = packages.speclab.orthopoly.jacobi;
opt = jac.defaults(varargin{:});

flags = false(0);
descriptions = {};

jint = jac.interval(opt);

% Gauss quadrature tests
[r,w] = jac.quad.gauss_quadrature(opt.N,opt);

descriptions{end+1} = 'Gauss quadrature nodes inside interval';
flags(end+1) = all(r<jint(2)) && all(r>jint(1));

descriptions{end+1} = 'Gauss quadrature satisfies polynomial integration accuracy';
% do scaled and shifted monomials -- otherwise roundoff errors might be large
flags(end+1) = true;
ps = jac.eval.eval_jacobi_poly(r,0:(2*opt.N-1),opt);
tol = 1e-10;
for q = 1:(2*opt.N-1)
  int_val = w.'*ps(:,q+1);
  flags(end) = and(flags(end),abs(int_val)<tol);
end

% Gauss-Radau nodes
opt.r = opt.scale+opt.shift;
[r,w] = jac.quad.gauss_radau_quadrature(opt.N,opt);

descriptions{end+1} = 'Gauss-Radau quadrature nodes inside interval';
flags(end+1) = all(r<=jint(2)+tol) && all(r>=jint(1)-tol);

descriptions{end+1} = 'Gauss Radau point at +1';
flags(end+1) = abs(r(end)-jint(2))<tol;

descriptions{end+1} = 'Gauss-Radau quadrature satisfies polynomial integration accuracy';
% do scaled and shifted monomials -- otherwise roundoff errors might be large
flags(end+1) = true;
ps = jac.eval.eval_jacobi_poly(r,0:(2*opt.N-2),opt);
tol = 1e-10;
for q = 1:(2*opt.N-2)
  int_val = w.'*ps(:,q+1);
  flags(end) = and(flags(end),abs(int_val)<tol);
end

opt.r = -opt.scale+opt.shift;
[r,w] = jac.quad.gauss_radau_quadrature(opt.N,opt);

descriptions{end+1} = 'Gauss-Radau quadrature nodes inside interval';
flags(end+1) = all(r<=jint(2)+tol) && all(r>=jint(1)-tol);

descriptions{end+1} = 'Gauss Radau point at -1';
flags(end+1) = abs(r(1)-jint(1))<tol;

descriptions{end+1} = 'Gauss-Radau quadrature satisfies polynomial integration accuracy';
% do scaled and shifted monomials -- otherwise roundoff errors might be large
flags(end+1) = true;
ps = jac.eval.eval_jacobi_poly(r,0:(2*opt.N-2),opt);
tol = 1e-10;
for q = 1:(2*opt.N-2)
  int_val = w.'*ps(:,q+1);
  flags(end) = and(flags(end),abs(int_val)<tol);
end

opt = rmfield(opt,'r');
opt.r1 = jint(1);
opt.r2 = jint(2);
[r,w] = jac.quad.gauss_lobatto_quadrature(opt.N,opt);

descriptions{end+1} = 'Gauss-Lobatto quadrature nodes inside interval';
flags(end+1) = all(r<=jint(2)+tol) && all(r>=jint(1)-tol);

descriptions{end+1} = 'Gauss Lobatto points at -1,+1';
flags(end+1) = (abs(r(1)-jint(1))<tol) && ...
               (abs(r(end)-jint(2))<tol);

descriptions{end+1} = 'Gauss-Lobatto quadrature satisfies polynomial integration accuracy';
% do scaled and shifted monomials -- otherwise roundoff errors might be large
flags(end+1) = true;
ps = jac.eval.eval_jacobi_poly(r,0:(2*opt.N-3),opt);
tol = 1e-10;
for q = 1:(2*opt.N-3)
  int_val = w.'*ps(:,q+1);
  flags(end) = and(flags(end),abs(int_val)<tol);
end

% Mass matrix: use gauss quadrature
[r,w] = jac.quad.gauss_quadrature(opt.N,opt);

ps = jac.eval.eval_jacobi_poly(r,0:(opt.N-1),opt);
mass = ps'*spdiags(w,0,opt.N,opt.N)*ps;

descriptions{end+1} = 'Gauss-quadrature mass matrix is diagonal';
flags(end+1) = norm(mass - opt.scale*eye(opt.N));

% Connection coefficient building blocks
r = linspace(jint(1),jint(2),1000).';
Nr = length(r);
promote_opt = opt;
promote_opt.alpha = promote_opt.alpha+1;
ns = 0:(promote_opt.N-1);
ns_1 = 0:(promote_opt.N);
ns_2 = 0:(promote_opt.N+1);

ps = jac.eval.eval_jacobi_poly(r,ns, promote_opt);
mu = jac.coefficients.one_minus_r_times_p(ns,promote_opt.alpha,promote_opt.beta,promote_opt);
ps_temp = jac.eval.eval_jacobi_poly(r,ns_2,opt);

tempopt = opt;
tempopt.beta = 0;
tempopt.alpha = 1;
err = spdiags(jac.weights.weight(r,tempopt),0,Nr,Nr)*ps - ps_temp(:,1:opt.N)*spdiags(mu(:,1),0,opt.N,opt.N) - ...
                                ps_temp(:,2:(opt.N+1))*spdiags(mu(:,2),0,opt.N,opt.N);

descriptions{end+1} = '(1-r)*P coefficients';
tol = 5e-3;
flags(end+1) = norm(err)<tol;

promote_opt = opt;
promote_opt.beta = promote_opt.beta+1;
ps = jac.eval.eval_jacobi_poly(r,ns, promote_opt);
mu = jac.coefficients.one_plus_r_times_p(ns,promote_opt.alpha,promote_opt.beta,promote_opt);
tempopt.alpha = 0;
tempopt.beta = 1;
err = spdiags(jac.weights.weight(r,tempopt),0,Nr,Nr)*ps - ps_temp(:,1:opt.N)*spdiags(mu(:,1),0,opt.N,opt.N) - ...
                                ps_temp(:,2:(opt.N+1))*spdiags(mu(:,2),0,opt.N,opt.N);

descriptions{end+1} = '(1+r)*P coefficients';
tol = 5e-3;
flags(end+1) = norm(err)<tol;

promote_opt = opt;
promote_opt.beta = promote_opt.beta+1;
promote_opt.alpha = promote_opt.alpha+1;
tempopt.alpha = 1;
tempopt.beta = 1;
ps = jac.eval.eval_jacobi_poly(r,ns, promote_opt);
mu = jac.coefficients.one_minus_r_squared_times_p(ns,promote_opt.alpha,promote_opt.beta,promote_opt);

err = spdiags(jac.weights.weight(r,tempopt),0,Nr,Nr)*ps - ps_temp(:,1:opt.N)*spdiags(mu(:,1),0,opt.N,opt.N) - ...
                                   ps_temp(:,2:(opt.N+1))*spdiags(mu(:,2),0,opt.N,opt.N) - ...
                                   ps_temp(:,3:end)*spdiags(mu(:,3),0,opt.N,opt.N); 

descriptions{end+1} = '(1-r^2)*P coefficients';
tol = 5e-3;
flags(end+1) = norm(err)<tol;

%%%%% Test modal stuff
% Test function, for any scale,shift:
f = @(x) sin(5*x);
df = @(x) 5*cos(5*x);

[r,w] = jac.quad.gauss_quadrature(opt.N,opt);
ps = jac.eval.eval_jacobi_poly(r,0:(opt.N-1), opt);
gauss_modes = ps'*spdiags(w,0,opt.N,opt.N)*f(r);

r_refined = linspace(jint(1),jint(2),1000).';
ps_refined = jac.eval.eval_jacobi_poly(r_refined,0:(opt.N-1),opt);
dopt = opt;
dopt.d = 1;
dps_refined = jac.eval.eval_jacobi_poly(r_refined,0:(dopt.N-1),dopt);

descriptions{end+1} = 'Gauss interpolation error';
tol = 1e-4;
flags(end+1) = max(abs(ps_refined*gauss_modes - f(r_refined)))<tol;

descriptions{end+1} = 'Gauss interpolation-differentiation error';
tol = 5e-2;
flags(end+1) = max(abs(dps_refined*gauss_modes - df(r_refined)))<tol;

A = ceil(opt.alpha);
B = ceil(opt.beta);
demoted_opt = opt;
demoted_opt.alpha = opt.alpha - A;
demoted_opt.beta = opt.beta - B;
[r_demoted,w_demoted] = jac.quad.gauss_quadrature(demoted_opt.N,demoted_opt);
ps_demoted = jac.eval.eval_jacobi_poly(r_demoted,0:(opt.N-1),demoted_opt);
modes = ps_demoted'*spdiags(w_demoted,0,opt.N,opt.N)*f(r_demoted);

C = jac.connection.integer_separation_connection_matrix(opt.N,demoted_opt.alpha,demoted_opt.beta,A,B);
promoted_modes = C*modes;

descriptions{end+1} = 'Integer connection matrix';
tol = 1e-5;
flags(end+1) = max(abs(gauss_modes(1:(end-A-B)) - ...
                       promoted_modes(1:(end-A-B))))<tol;

varargout{1} = flags;
varargout{2} = descriptions;
