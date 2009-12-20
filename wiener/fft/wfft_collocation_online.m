function[F] = wfft_collocation_online(f,data)
% [F] = wfft_collocation_online(f,data)
%
%     Computes the Wiener collocation FFT of the function evaluations f using
%     the precomputed data stored in the struct data.

persistent ffft_online
if isempty(ffft_online)
  from speclab.fourier.fft import ffft_online
end

f = f./data.weight;
F = ffft_online(f,data.fourier_data);
