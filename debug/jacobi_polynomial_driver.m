function[varargout] = jacobi_polynomial_debug(varargin)
% [FLAGS,DESCRIPTIONS] = JACOBI_POLYNOMIAL_DRIVER({N=10,ALPHA=-1/2,BETA=-1/2,SCALE=1,SHIFT=0})
% 
%     Runs various tests on speclab's Jacobi polynomial package. FLAGS is a
%     boolean array containing success indicators for each test, whose
%     description is in the cell array DESCRIPTIONS.

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
opt = jac.defaults(varargin{:});

flags = [];
descriptions = {};

jint = jac.interval(opt);

% Gauss quadrature tests
[r,w] = jac.quad.gauss_quadrature(opt.N,opt);

descriptions{end+1} = "Gauss quadrature nodes inside interval";
flags(end+1) = all(r<interval(2)) && all(r>interval(1));

descriptions{end+1} = "Gauss quadrature satisfies polynomial integration accuracy";
% do scaled and shifted monomials -- otherwise roundoff errors might be large
flags(end+1) = true;
tol = 1e-12;
for q = 1:(2*opt.N-1)
  int_val = w.'*((r-opt.shift)/opt.scale).^q;
  flags(end) = and(flags(end),abs(int_val)<tol);
end

% Gauss-Radau nodes
opt.r = 1;
[r,w] = jac.quad.gauss_radau_quadrature(opt.N,opt);

descriptions{end+1} = "Gauss-Radau quadrature nodes inside interval";
flags(end+1) = all(r<=interval(2)) && all(r>=interval(1));

descriptions{end+1} = "Gauss Radau point at +1";
flags(end+1) = abs(r(end)-interval(2))<tol;

descriptions{end+1} = "Gauss-Radau quadrature satisfies polynomial integration accuracy";
% do scaled and shifted monomials -- otherwise roundoff errors might be large
flags(end+1) = true;
tol = 1e-12;
for q = 1:(2*opt.N-2)
  int_val = w.'*((r-opt.shift)/opt.scale).^q;
  flags(end) = and(flags(end),abs(int_val)<tol);
end

opt.r = -1;
[r,w] = jac.quad.gauss_radau_quadrature(opt.N,opt);

descriptions{end+1} = "Gauss-Radau quadrature nodes inside interval";
flags(end+1) = all(r<=interval(2)) && all(r>=interval(1));

descriptions{end+1} = "Gauss Radau point at -1";
flags(end+1) = abs(r(1)-interval(1))<tol;

descriptions{end+1} = "Gauss-Radau quadrature satisfies polynomial integration accuracy";
% do scaled and shifted monomials -- otherwise roundoff errors might be large
flags(end+1) = true;
tol = 1e-12;
for q = 1:(2*opt.N-2)
  int_val = w.'*((r-opt.shift)/opt.scale).^q;
  flags(end) = and(flags(end),abs(int_val)<tol);
end

opt = rmfield(opt,'r');
[r,w] = jac.quad.gauss_lobatto_quadrature(opt.N,opt);

descriptions{end+1} = "Gauss-Lobatto quadrature nodes inside interval";
flags(end+1) = all(r<=interval(2)) && all(r>=interval(1));

descriptions{end+1} = "Gauss Lobatto points at -1,+1";
flags(end+1) = (abs(r(1)-interval(1))<tol) && ...
               (abs(r(end)-interval(2))<tol);

descriptions{end+1} = "Gauss-Lobatto quadrature satisfies polynomial integration accuracy";
% do scaled and shifted monomials -- otherwise roundoff errors might be large
flags(end+1) = true;
tol = 1e-12;
for q = 1:(2*opt.N-3)
  int_val = w.'*((r-opt.shift)/opt.scale).^q;
  flags(end) = and(flags(end),abs(int_val)<tol);
end

% Mass matrix: use gauss quadrature
[r,w] = jac.quad.gauss_quadrature(opt.N,opt);

ps = jac.eval.eval_jacobi_poly(r,0:(opt.N-1),opt);
mass = ps'*spdiags(w,0,opt.N,opt.N)*ps;

descriptions{end+1} = "Gauss-quadrature mass matrix is diagonal";
flags(end+1) = norm(mass - opt.scale*eye(opt.N));
