function[f] = iwfft_galerkin_online(F,data)
% [f] = iwfft_galerkin_online(F,data)
%
%     Computes the iwfft 'Galerkin' algorithm using an online/offline
%     decomposition. The input data comes from wfft_galerkin_overhead.

persistent wiener fourier wconnect
if isempty(wiener)
  from speclab import wiener
  from speclab import fourier
  from speclab.fourier import negative_integer_separation_connection_online as wconnect
end
%global packages;
%wiener = packages.speclab.wiener;
%fourier = packages.speclab.fourier;
%wconnect = fourier.connection.negative_integer_separation_connection_online;

if data.S>0
  F = wconnect(F,data.connection);
end
F = wiener.operators.wiener_weight_multiply(F,data.opt);
f = fourier.fft.iffft_online(F,data.fftdata);
