function[F] = jfft_online(f,data)
% [F] = jfft_online(f,data)
%
%     Uses precomputed data from jfft_overhead to perform the Jacobi FFT.

persistent chebfft_online
if isempty(chebfft_online)
  from speclab.orthopoly1d.jacobi.fft import chebfft_online
end
%global packages;
%jac = packages.speclab.orthopoly1d.jacobi;

F = chebfft_online(f,data.chebdata);
F = data.C*F;
