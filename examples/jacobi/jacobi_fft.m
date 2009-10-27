%%% Example script for Jacobi approximation
%%% Jacobi FFT: usage of the FFT for Jacobi polynomial expansions

% If 2*alpha and 2*beta are odd integers, the FFT may be used to compute Jacobi
% polynomial expansion modal coefficients in N*log(N) time. 
clear
jac = from_package_import_as('speclab.orthopoly1d', 'jacobi');
%jac = packages.speclab.orthopoly1d.jacobi;
M = 100;  % # of times to repeat calculation for timing

% Unfortunately, the FFT is not supported if the conditions on alpha and beta
% are not met.
jopt.alpha = 3;
jopt.beta = 2;
copt.alpha = -1/2;  % Since -1/2 is default, this doesn't need to be set explicity
copt.beta = -1/2;
Nq = 800; % Quadrature size
N = 800;  % # of modal coefficients
[cr,cw] = jac.quad.gauss_quadrature(Nq,copt); % Need the Chebyshev nodes
fr = zeros(size(cr));
fprintf('Try to do FFT with non-conforming alpha,beta: alpha = %2.1f, beta = %2.1f \n', jopt.alpha, jopt.beta);
try
  jac.fft.jfft(fr,jopt);
catch
  % The following error was given to Matlab:
  fprintf('Cannot use the FFT when 2*alpha and 2*beta are not odd integers\n');
end

% Now for a real example
jopt.alpha = 11/2;
jopt.beta = 7/2;
fprintf('\nNow using alpha = %2.1f, beta = %2.1f\n\n', jopt.alpha, jopt.beta);
% You don't need to run the following for other applications: it's just a helper
% for this script. It determines (a) if a basis set can use the FFT, and (b)
% (alpha+1/2) and (beta+1/2)
[blah,A,B] = jac.fft.fftable(jopt);

% example function
f = @(r) atan(50*r);

% First let's use direct Gauss quadrature: takes a *long* time to compute overhead
% Most of the computational time required in this script is in the following
% tic/toc-surrounded lines:
tic;
  [jr,jw] = jac.quad.gauss_quadrature(Nq,jopt);
  jr_vandermonde = jac.eval.eval_jacobi_poly(jr,0:(N-1),jopt);
  jr_vandermonde_inverse = jr_vandermonde'*spdiags(jw,0,Nq,Nq);
gaussquad_overhead_time = toc;
fjr = f(jr);
% Direct quadrature modes:
guassquad_modes = zeros([N 1]);
tic;
  for q=1:M
    gaussquad_modes = jr_vandermonde_inverse*fjr;
  end
gaussquad_online_time = toc/M;

% Now let's use the Chebyshev rule:
tic;
  [cr,cw] = jac.quad.gauss_quadrature(Nq,copt);
  cr_vandermonde = jac.eval.eval_jacobi_poly(cr,0:(N-1),jopt);
  cw = cw.*jac.weights.weight(cr,'alpha',A,'beta',B);
  cr_vandermonde_inverse = cr_vandermonde'*spdiags(cw,0,Nq,Nq);
chebquad_overhead_time = toc;
fcr = f(cr);
chebquad_modes = zeros([N 1]);
tic;
  for q=1:M
    chebquad_modes = cr_vandermonde_inverse*fcr;
  end
chebquad_online_time = toc/M;
% chebquad_time should be about the same as gaussquad_time

% Let's do the Jacobi FFT now: very little overhead is needed.
jfft_modes = zeros([N 1]);
tic;
for q = 1:M
  jfft_modes = jac.fft.jfft(fcr,jopt);
end
jfft_all_time = toc/M;
% The above is slower than the previous two because jfft has to compute O(N)
% overhead that is written in an m-file: much slower than precompiled BLAS 2
% operations that is matlab's * operator.

% Let's do the FFT with an online/offline decomposition
fftdata = jac.fft.jfft_overhead(Nq,jopt);
tic;
for q = 1:M
  jfft_modes = jac.fft.jfft_online(fcr,fftdata);
end
jfft_online_time = toc/M;

% You may check yourself that chebquad_modes, gaussquad_modes, and jfft_modes
% are nearly identical. For large alpha, beta, using the gaussquad_modes method
% is inaccurate due to roundoff errors in evaluating the quadrature rule for
% Jacobi polynomials

% Now output information
fprintf('***  Direct Gauss quadrature timing  ***\n')
fprintf('            Overhead time: %1.4e\n', gaussquad_overhead_time);
fprintf('            Online   time: %1.4e\n', gaussquad_online_time);

fprintf('***  Chebyshev-Gauss quadrature timing  ***\n')
fprintf('            Overhead time: %1.4e\n', chebquad_overhead_time);
fprintf('            Online   time: %1.4e\n', chebquad_online_time);

fprintf('***  Jacobi FFT timing  ***\n')
fprintf('            Overhead time: %1.4e\n', jfft_all_time - jfft_online_time);
fprintf('            Online   time: %1.4e\n', jfft_online_time);

% Note: the Jacobi FFT is in theory *much* faster than this example shows: the
% quadrature timings are all for Matlab's precompiled BLAS operations. The FFT
% times are dominated by m-file runtime-compiled code.

% The "ijfft" routines are used in an identical way. However, the fftdata given
% as input to ijfft_online is identical to that given to jfft_online. 
