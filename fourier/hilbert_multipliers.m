function[F] = hilbert_multipliers(N)
% hilbert_multipliers -- Frequency-space multipliers for the (periodic) Hilbert transform
%
% F = hilbert_multipliers(N)
%
%      For an N-mode Fourier expansion, this function returns the length-N
%      vector of Fourier-space multipliers that correspond to taking the Hilbert
%      transform.

persistent integer_range
if isempty(integer_range)
  from speclab.common import integer_range
end

ks = integer_range(N);
F = zeros(size(ks));

F(ks>0) = -i;
F(ks<0) = i;
