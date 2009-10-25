%%% Example script for Wiener approximation
%%% Stiffness matrix: the Wiener stiffness matrix is sparse, O(N) spectral
%%% radius

clear
global packages;
wiener = packages.speclab.wiener;
ltex = packages.common.typelatex;
explot = packages.common.explot;

opt.s = 1;
N = 130;
stiffmat_s1 = wiener.matrices.wiener_stiffness_matrix(N,opt);

opt.s = -1/2 + 10*rand;
stiffmat_s = wiener.matrices.wiener_stiffness_matrix(N,opt);

subplot(1,2,1);
spy(stiffmat_s1);
axis('square');
ltex(title('Sparsity pattern, $s=1$'));

subplot(1,2,2);
spy(stiffmat_s);
axis('square');
ltex(title(['Sparsity pattern, $s=' num2str(opt.s) '$']));

% Calculate spectral radius vs. N
Ns = N:-2:2; % N must be even for this for-loop to work
s1_radius = zeros(size(Ns));
s1_temp = stiffmat_s1;
s_radius = zeros(size(Ns));
s_temp = stiffmat_s;
for n=Ns
  s1_radius(n/2) = max(abs(eig(full(s1_temp))));
  s_radius(n/2) = max(abs(eig(full(s_temp))));

  s1_temp = s1_temp(2:end-1,2:end-1);
  s_temp = s_temp(2:end-1,2:end-1);
end
s1_radius = fliplr(s1_radius);
s_radius = fliplr(s_radius);

figure;
explot(Ns,s1_radius, 'k', Ns, s_radius, 'r--');
ltex(legend('$s=1$', ['$s=' num2str(opt.s) '$']));
ltex(xlabel('$N$'));
ltex(ylabel('$\rho(S)$'));
ltex(title('Spectral radius of $N \times N$ stiffness matrix'));

% The above really is the stiffness matrix:
[x,w] = wiener.quad.pi_gauss_quadrature(2*N);
ks = packages.speclab.common.integer_range(N);
ws = wiener.eval.wiener_function(x,ks);
dws = wiener.eval.derivative_wiener_function(x,ks);

stiffmat_s1_numerical = ws'*spdiags(w,0,2*N,2*N)*dws;

% Unless s is an integer, the s=1 quadrature rule isn't that accurate
[x,w] = wiener.quad.pi_gauss_quadrature(2*N,opt);
ws = wiener.eval.wiener_function(x,ks,opt);
dws = wiener.eval.derivative_wiener_function(x,ks,opt);

stiffmat_s_numerical = ws'*spdiags(w,0,2*N,2*N)*dws;

err = norm(stiffmat_s1 - stiffmat_s1_numerical);
fprintf('s=1 stiffness matrix error is %1.4e\n', err);

err = norm(stiffmat_s - stiffmat_s_numerical);
fprintf('s=%1.3f stiffness matrix error is %1.4e\n', opt.s, err);

% You can check that the stiffness matrices are skew-Hermitian.
