%%% Example script for Jacobi approximation
%%% Gibbs resolution: re-expansion in Gegenbauer polynomials to remove Gibbs'
%%%                   oscillations

% This is a bit of a cheat: instead of working with the exact global modal
% coefficients, we work with quadrature on a crapload of local Chebyshev nodes.
% For this example, the difference is almost certainly negligible.
% You may want to take a look at the examples/fourier/intro.m file before going
% through this.

clear
from_package_import('speclab', 'fourier');
jac = from_package_import_as('speclab.orthopoly1d', 'jacobi');
[irange, ltex] = from_package_import_as('speclab.labtools', 'integer_range', 'typelatex')
ltex = 
%fourier = packages.speclab.fourier;
%jac = packages.speclab.orthopoly1d.jacobi;
%irange = packages.speclab.common.integer_range;
%ltex = packages.labtools.typelatex;

% Let's try to approximate a discontinuous function by a Fourier Series
f = @(x) (x>=0).*(x<=1).*exp(-x) + ...
         (x>1).*(x<=2).*((x/2).^8) + ...
         (x>2).*(x<=3).*(sin(3*x)) + ...
         (x>3).*(x<=4).*(cos(sqrt(x-3)));

Nq = 1e5;  % Fourier quadrature is cheap with the FFT
N = 250;   % Let's not take so many modes
fopt = fourier.affine_scaling([0,4]);
[theta,w] = fourier.quad.gauss_quadrature(Nq,fopt);
[theta_N,w] = fourier.quad.gauss_quadrature(N,fopt);

fourier_modes = fourier.fft.ffft(f(theta_N),fopt);
% Meh, this is techical...I'm just zero-padding the modal coefficient vector to
% have Nq entries
fourier_modes_padded = zeros([Nq 1]);
ks = irange(N) + 1 - min(irange(Nq));
fourier_modes_padded(ks) = fourier_modes;
kmin = min(irange(N));

% Now fourier_interp has crazy Gibbs oscillations
fourier_interp = fourier.fft.iffft(fourier_modes_padded,fopt);

subplot(2,2,1); 
plot(theta, f(theta));
title('Original function');
ltex(xlabel('$\theta$'));
axis([0,4,-0.5 1.25]);
subplot(2,2,2);
plot(theta, real(fourier_interp)); % See plot: ugly-looking, eh?
title('Gibbs'' oscillation-polluted Fourier approximation');
ltex(xlabel('$\theta$'));
axis([0,4,-0.5 1.25]);

beta_factor = 0.10;  % See Gottlieb/Shu for explanations of what these are
alpha_factor = 0.18;
gegenbauer_postprocessing = zeros(size(fourier_interp));

% Assume we know the sub-intervals where the function is analytic:
intervals = [0, 1;...
             1, 2;...
             2, 3;...
             3, 4];

Nplots = 20:20:200;  % For these #'s of Fourier modes, we'll do the method
Linf_error = zeros(size(Nplots));
NNcount = 1;
Nq = round(Nq/100);  %  Don't really need 10^5 points now, do we?
for NN=Nplots
  ks = irange(NN);
  NN_fourier_modes = fourier_modes(ks+1-kmin);
  % m and lambda defined by Gottlieb & Shu
  lambda = round(alpha_factor*floor(NN/2));
  m = round(beta_factor*floor(NN/2));

  for q = 1:size(intervals,1) % for each sub-interval
    interval = intervals(q,:);
    copt = jac.affine_scaling(interval);
    % Local chebyshev nodes:
    [r,w] = jac.quad.gauss_quadrature(Nq,copt);

    % Interpolate modes to local Chebyshev modes: can't use the fft: use a
    % vandermonde matrix. This line looping over q, NN does the same thing lots
    % of times and can be optimized.
    vandermonde = fourier.eval.fseries(r,ks,fopt);
    % This is the Gibbs'-polluted solution evaluated at the Chebyshev points:
    fr = real(vandermonde*NN_fourier_modes);

    % Use the Jacobi fft to get Gegenbauer modes
    jopt = copt; % for scale/shift factors
    jopt.alpha = -1/2 + lambda;
    jopt.beta = jopt.alpha;

    gegenbauer_modes = jac.fft.jfft(fr,jopt);
    gegenbauer_modes = gegenbauer_modes(1:m); % Only take the first m modes

    % find theta's inside interval
    theta_inds = interval(1)<theta & theta<interval(2);
    temp = theta(theta_inds);
    % Interpolate Gegenbauer approximation of order beta_factor*NN to these
    % locations--again, no FFT possible.
    vandermonde = jac.eval.eval_jacobi_poly(temp,0:(m-1),jopt);
    gegenbauer_postprocessing(theta_inds) = vandermonde*gegenbauer_modes;
  end

  % What's the Linf error for this value of NN?
  Linf_error(NNcount) = max(abs(gegenbauer_postprocessing - f(theta)));
  NNcount = NNcount + 1;
end

subplot(2,2,3); 
plot(theta, gegenbauer_postprocessing); 
ltex(xlabel('$\theta$'));
title('Gegenbauer post-processed solution for N=200');
axis([0,4,-0.5 1.25]);
subplot(2,2,4);
semilogy(Nplots, Linf_error, 'b.-'); 
ltex(xlabel('$N$: Number of Fourier modes'));
ltex(ylabel('$L^\infty$ error'));
ltex(title('$L^\infty$ convergence of Gegenbauer postprocessing'));
