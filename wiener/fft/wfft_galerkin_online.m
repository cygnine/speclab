function[F] = wfft_galerkin_online(f,data)
% [F] = wfft_galerkin_online(f,data)
%
%     Computes the wfft 'Galerkin' algorithm using an online/offline
%     decomposition. The input data comes from wfft_galerkin_overhead.

persistent wiener fourier wconnect
if isempty(wiener)
  from speclab import wiener fourier
  from speclab.fourier.connection import positive_integer_separation_connection_online as wconnect
end

F = fourier.fft.ffft_online(f,data.fftdata);
F = wiener.operators.wiener_weight_divide(F,data.opt);
if data.S>0
  F = wconnect(F,data.connection);
  F(end-data.ST:end) = 0;
end

F(1:(data.ST+1)) = 0;
