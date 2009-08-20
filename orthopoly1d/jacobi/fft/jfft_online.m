function[F] = jfft_online(f,data)
% [F] = jfft_online(f,data)
%
%     Uses precomputed data from jfft_overhead to perform the Jacobi FFT.

global handles;
jac = handles.speclab.orthopoly1d.jacobi;

F = jac.fft.chebfft_online(f,data.chebdata);
F = data.C*F;
