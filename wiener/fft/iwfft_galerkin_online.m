function[f] = iwfft_galerkin_online(F,data)
% [f] = iwfft_galerkin_online(F,data)
%
%     Computes the iwfft 'Galerkin' algorithm using an online/offline
%     decomposition. The input data comes from wfft_galerkin_overhead.

persistent wiener_weight_multiply iffft_online wconnect
if isempty(wconnect)
  from speclab.wiener.operators import wiener_weight_multiply
  from speclab.fourier.fft import iffft_online
  from speclab.fourier.connection import negative_integer_separation_connection_online as wconnect
end

if data.S>0
  F = wconnect(F,data.connection);
end
F = wiener_weight_multiply(F,data.opt);
f = iffft_online(F,data.fftdata);
