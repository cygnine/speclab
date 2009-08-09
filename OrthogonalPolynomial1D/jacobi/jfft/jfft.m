function[F] = jfft(f,varargin)
% [F] = jfft(f,{points='gq', alpha=-1/2, beta=-1/2, normalization='normal', scale=1})
%
%     If 2*alpha and 2*beta are odd, this function computes the Jacobi
%     Polynomial modal coefficients of the function point-values f. The
%     locations of the nodes must be a **Chebyshev** fft-compatible set of
%     points. See chebfft. 

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;

tol = 1e-12;
A = opt.alpha + 1/2;
B = opt.beta + 1/2;

if abs(A-round(A))>tol || abs(B - round(B))>tol
  error('This basis type is not fftable');
end

A = round(A);
B = round(B);
N = size(f,1);

F = jac.jfft.chebfft(f,opt);
C = jac.connection.integer_separation_connection_matrix(N,opt.alpha,opt.beta,A,B);
F = C*F;
