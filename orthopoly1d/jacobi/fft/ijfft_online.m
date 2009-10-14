function[f] = jifft_online(F,data)
% [f] = jifft_online(F,data)
%
%     Uses precomputed data from jifft_overhead to perform the Jacobi IFFT.

global handles;
jac = handles.speclab.orthopoly1d.jacobi;
la = handles.common.linalg;

F = la.triu_sparse_invert(data.C,F,'bandwidth',data.A+data.B+1);
f = jac.fft.chebifft_online(F,data.chebdata);