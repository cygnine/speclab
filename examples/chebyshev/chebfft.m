%%% Example script for Chebyshev approximation
%%% fft: how to use the FFT for nodal/modal transformations

% Computing modal coefficients in a Chebyshev expansion can be done using the
% FFT if the nodal evaluations are placed at Chebyshev-Gauss-like points. At
% present, only Chebyshev-Gauss is implemented. Radau/Lobatto is TODO. 
% Matlab is so slow at user-programmed routines that using the FFT is really
% only necessary if N is too large to invert/store a Vandermonde matrix (~2000).
% Otherwise, one might as well store the inverted Vandermonde matrix and apply
% it as regular matrix-multiplication; it'll be much faster than the FFT. 

clear
from speclab.orthopoly1d import jacobi as cheb
%cheb = packages.speclab.orthopoly1d.jacobi;
map = cheb.affine_scaling([3,4]); % for no particularly compelling reason

N = 700;
M = 100; % # of times to run FFT test
[r,w] = cheb.quad.gauss_quadrature(N,map);
vandermonde = cheb.eval.eval_jacobi_poly(r,0:(N-1),map);

f = @(r) exp(sin(20*r));
fr = f(r);

% The first run is a wash because matlab has to store it's own overhead, so
% carry out each task M times and take average
modal_coefficients_1 = zeros([N,1]);
modal_coefficients_2 = zeros([N,1]);
modal_coefficients_3 = zeros([N,1]);
vinv = vandermonde'*spdiags(w,0,N,N);

vinv*fr;
tic;
for q = 1:M;
  modal_coefficients_1 = vinv*fr;
end
time1 = toc;

cheb.fft.chebfft(fr);
tic;
for q = 1:M
  modal_coefficients_2 = cheb.fft.chebfft(fr);
end
time2 = toc;

err = norm(modal_coefficients_1 - modal_coefficients_2);
fprintf('FFT error is %1.3e\n', err);

% cheb_fft requires a bit of overhead since the Chebyshev polys aren't exactly
% direct maps of the canonical Fourier basis. You can store this overhead for
% later use:
fftdata = cheb.fft.chebfft_overhead(N,map);
cheb.fft.chebfft_online(fr,fftdata);
tic;
for q = 1:M;
  modal_coefficients_3 = cheb.fft.chebfft_online(fr,fftdata);
end
time3 = toc;

str = ['Direct quadrature timing:   %1.3e\n' ...
       'Direct FFT timing:          %1.3e\n' ...
       'Overhead/online FFT timing: %1.3e\n'];
fprintf(str, time1/M, time2/M, time3/M);
