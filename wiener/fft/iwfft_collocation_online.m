function[f] = iwfft_collocation_online(F,data)
% [f] = iwfft_collocation_online(F,data)
%
%     Computes the Wiener collocation IFFT of the modal coefficients F using the
%     precomputed data stored in the struct data.

global packages;
fourier = packages.speclab.fourier;

f = fourier.fft.iffft_online(F,data.fourier_data);
f = f.*data.weight;
