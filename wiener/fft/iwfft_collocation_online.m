function[f] = iwfft_collocation_online(F,data)
% [f] = iwfft_collocation_online(F,data)
%
%     Computes the Wiener collocation IFFT of the modal coefficients F using the
%     precomputed data stored in the struct data.

persistent iffft_online
if isempty(iffft_online)
  from speclab.fourier.fft import iffft_online
end

f = iffft_online(F,data.fourier_data);
f = f.*data.weight;
