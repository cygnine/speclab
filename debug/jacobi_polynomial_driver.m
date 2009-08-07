function[varargout] = jacobi_polynomial_driver(varargin)
% [FLAGS,DESCRIPTIONS] = JACOBI_POLYNOMIAL_DRIVER({N=10,ALPHA=-1/2,BETA=-1/2,SCALE=1,SHIFT=0})
% 
%     Runs various tests on speclab's Jacobi polynomial package. FLAGS is a
%     boolean array containing success indicators for each test, whose
%     description is in the cell array DESCRIPTIONS.

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
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
tol = 1e-12;
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
tol = 1e-12;
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
tol = 1e-12;
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
tol = 1e-12;
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

varargout{1} = flags;
varargout{2} = descriptions;
