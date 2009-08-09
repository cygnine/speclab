function[f] = jifft(F,varargin)
% [f] = jifft(F,{points='gq', alpha=-1/2, beta=-1/2, normalization='normal', scale=1})
%
%     If 2*alpha and 2*beta are odd, this function computes the nodal
%     evaluations from the Jacobi Polynomial modal coefficients F using the
%     IFFT. The locations of the nodes must be a **Chebyshev** fft-compatible
%     set of points. See chebifft. 

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
la = handles.common.linalg;

tol = 1e-12;
A = opt.alpha + 1/2;
B = opt.beta + 1/2;

if abs(A-round(A))>tol || abs(B - round(B))>tol
  error('This basis type is not fftable');
end

A = round(A);
B = round(B);
N = size(f,1);

C = jac.connection.integer_separation_connection_matrix(N,opt.alpha,opt.beta,A,B);
F = la.triu_sparse_invert(C,F,'bandwidth',A+B+1);
f = jac.jfft.chebifft(F,opt);
