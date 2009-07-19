% Testing bare-bones stuff of Jacobi polynomials

clear;
global handles;
pss = handles.speclab.common.physical_scaleshift_1d;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;

N = ceil(100*rand(1));
alpha = -1/2 + 10*rand(1);
beta = -1/2 + 10*rand(1);
scale = 3*rand(1);
shift = randn(1);
r = -1;
r1 = -1;
r2 = 1;

% Test function:
tempf_shift = randn(1);
f = @(x) exp(-(x-tempf_shift).^2);
df = @(x) -2*(x-tempf_shift).*exp(-(x-tempf_shift).^2);

options = struct;
options.alpha = alpha;
options.beta = beta;
options.scale = scale;
options.shift = shift;
options.r = pss(r,options);
options.r1 = pss(r1,options);
options.r2 = pss(r2,options);
%options.normalization = 'normal';
options = jac.defaults(options);

chebcase = struct;
chebcase.scale = scale;
chebcase.shift = shift;

% Test basic evaluations
[x,w] = jac.gauss_quadrature(N,options);
[xfine,wfine] = jac.gauss_quadrature(10*N,chebcase);
n = 0:(N-1);
v = jac.eval_jacobi_poly(x,n,options);
vinv = v'*spdiags(w,0,N,N);
mass = vinv*v;

fprintf('##### Gauss quadrature ######\n');
fprintf('Error in mass matrix: %1.5e\n', norm(mass-diag(ones([N,1]))));

modes = vinv*f(x);
frec = jac.eval_jacobi_poly(xfine,n,options)*modes;
fprintf('Error in function reconstruction for %3d modes is %1.5e\n', N, norm(frec-f(xfine)));

options.d = 1;
% Test derivatives
dfrec = jac.eval_jacobi_poly(xfine,n,options)*modes;
fprintf('Error in derivative reconstruction for %3d modes is %1.5e\n\n', ...
  N, norm(wfine.*(dfrec-df(xfine))));

% Test the same junk with Gauss-Radau/Lobatto
fprintf('##### Gauss-Radau quadrature ######\n');
options.d = 0;
[x,w] = jac.gauss_radau_quadrature(N,options);
v = jac.eval_jacobi_poly(x,n,options);
vinv = v'*spdiags(w,0,N,N);
modes = vinv*f(x);
frec = jac.eval_jacobi_poly(xfine,n,options)*modes;
fprintf('Error in function reconstruction for %3d modes is %1.5e\n', N, norm(frec-f(xfine)));
options.d = 1;
dfrec = jac.eval_jacobi_poly(xfine,n,options)*modes;
fprintf('Error in derivative reconstruction for %3d modes is %1.5e\n\n', ...
  N, norm(wfine.*(dfrec-df(xfine))));

fprintf('##### Gauss-Lobatto quadrature ######\n');
options.d = 0;
[x,w] = jac.gauss_lobatto_quadrature(N,options);
v = jac.eval_jacobi_poly(x,n,options);
vinv = v'*spdiags(w,0,N,N);
modes = vinv*f(x);
frec = jac.eval_jacobi_poly(xfine,n,options)*modes;
fprintf('Error in function reconstruction for %3d modes is %1.5e\n', N, norm(frec-f(xfine)));
options.d = 1;
dfrec = jac.eval_jacobi_poly(xfine,n,options)*modes;
fprintf('Error in derivative reconstruction for %3d modes is %1.5e\n\n', ...
  N, norm(wfine.*(dfrec-df(xfine))));
