function[F] = hilbert_transform(f)
% hilbert_transform -- Computes the Hilbert transform via the FFT
%
% F = hilbert_transform(f)
%
%     Computes the Hilbert transform of f using the FFT.

persistent hilbert_multipliers ffft iffft spdiag
if isempty(ffft)
  from speclab.fourier.fft import ffft iffft
  from speclab.fourier import hilbert_multipliers
  from labtools import spdiag
end

N = size(f,1);

H = spdiag(hilbert_multipliers(N));

F = iffft(H*ffft(f));
