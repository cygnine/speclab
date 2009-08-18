function[F] = wfft_collocation_online(f,data)
% [F] = wfft_collocation_online(f,data)
%
%     Computes the Wiener collocation FFT of the function evaluations f using
%     the precomputed data stored in the struct data.

global handles;
fourier = handles.speclab.fourier;

f = f./data.weight;
F = fourier.fft.ffft_online(f,data.fourier_data);
