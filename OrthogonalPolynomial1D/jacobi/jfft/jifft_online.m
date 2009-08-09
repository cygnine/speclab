function[F] = jifft_online(f,data)
% [F] = jifft_online(f,data)
%
%     Uses precomputed data from jifft_overhead to perform the Jacobi IFFT.

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
la = handles.common.linalg;

F = la.triu_sparse_invert(C,F,'bandwidth',A+B+1);
f = jac.jfft.chebifft_online(F,data.chebdata);
