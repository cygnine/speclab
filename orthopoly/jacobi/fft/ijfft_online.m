function[f] = jifft_online(F,data)
% [f] = jifft_online(F,data)
%
%     Uses precomputed data from jifft_overhead to perform the Jacobi IFFT.

persistent chebifft_online triu_sparse_invert
if isempty(chebifft_online)
  from speclab.orthopoly.jacobi.fft import chebifft_online
  from labtools.linalg import triu_sparse_invert
end

F = triu_sparse_invert(data.C,F,'bandwidth',data.A+data.B+1);
f = chebifft_online(F,data.chebdata);
