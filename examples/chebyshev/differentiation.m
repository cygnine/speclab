%%% Example script for Chebyshev approximation
%%% Differentiation: how to take derivatives

clear
from speclab.orthopoly import jacobi as cheb
%cheb = packages.speclab.orthopoly.jacobi;

map = cheb.affine_scaling([-pi/2,exp(1)]);  % why not?

f0 = @(r) exp(sin(r));
f1 = @(r) cos(r).*exp(sin(r));
f2 = @(r) -sin(r).*exp(sin(r)) + cos(r).^2.*exp(sin(r));
N = 200;

% Let's do Gauss-Radau quadrature
% Radau quadrature takes an optional input 'r' that specifies where the fixed
% quadrature node goes. Since it's optional, we'll bundle it in with the rest of
% the optional arguments.
map.r = -pi/2;
[r,w] = cheb.quad.gauss_radau_quadrature(N,map);

% Giving the (meaningless) optional input map.r to eval_jacobi_poly doesn't
% matter: it'll just ignore it.
vandermonde = cheb.eval.eval_jacobi_poly(r,0:(N-1),map);

% To evaluate derivatives in orthogonal polynomial evaluation routines, the
% optional input d needs to be specified.
map.d = 1;
first_derivatives = cheb.eval.eval_jacobi_poly(r,0:(N-1),map);
map.d = 2;
second_derivatives = cheb.eval.eval_jacobi_poly(r,0:(N-1),map);
% You can set d to any whole number (0,1,2,...), but I wouldn't numerically
% trust anything bigger than 4.

% Explicitly compute a differentiation matrix
vandermonde_inverse = vandermonde'*spdiags(w,0,N,N); 
% (the inverse is only approximate since it's less-accurate Radau quadrature)
diff1mat = first_derivatives*vandermonde_inverse;
diff2mat = second_derivatives*vandermonde_inverse;

fr = f0(r);
err = norm(diff1mat*fr - f1(r));
fprintf('Using direct quadrature:\n');
fprintf('First derivative error is %1.3e\n', err);
err = norm(diff2mat*fr - f2(r));
fprintf('Second derivative error is %1.3e\n', err);

%% O(N log N) differentiation
% Using the FFT, we can present an asymptotically O(N log N) differentiation
% method. 
% 1.) Compute modes - FFT O(N log N)
% 2.) `Promote' Jacobi (-1/2,-1/2) modes to (1/2,1/2) with the derivative - O(N)
% 3.) `Demote' Jacobi (1/2,1/2) modes back down to (-1/2, -1/2) by three-term
%     recurrence relation/inversion of tridiagonal system - O(N)
% 4.) Compute nodal evaluations - IFFT O(N log N)
%
% Since the chebfft is supported only on the Gauss nodes (at present), we have
% to go back to Gauss nodes
[r,w] = cheb.quad.gauss_quadrature(N,map);
fr = f0(r);

chebmodes = cheb.fft.chebfft(fr,map);
diff1modes = cheb.operators.stiffness_operator(chebmodes,map);
diff2modes = cheb.operators.stiffness_operator(diff1modes,map);

err = norm(cheb.fft.chebifft(diff1modes) - f1(r));
fprintf('\nUsing the FFT + stiffness matrix application:\n');
fprintf('First derivative error is %1.3e\n', err);
err = norm(cheb.fft.chebifft(diff2modes) - f2(r));
fprintf('Second derivative error is %1.3e\n', err);

% One can see that using the latter method yields more accurate results in this
% case. This difference is *not* apparent due to the more accurate Gauss rule in
% the second case.
