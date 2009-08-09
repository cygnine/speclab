function[data] = jfft_overhead(N,varargin)
% [data] = jfft_overhead(N,{points='gq', alpha=-1/2, beta=-1/2, normalization='normal', scale=1})
%
%     If 2*alpha and 2*beta are odd, this function computes the overhead data
%     necessary for using the N-point Jacobi FFT.

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
inputs = {'points', 'alpha', 'beta', 'normalization', 'scale'};
defaults = {'gq', -1/2, -1/2, 'normal', 1};

tol = 1e-12;
A = opt.alpha + 1/2;
B = opt.beta + 1/2;

if abs(A-round(A))>tol || abs(B - round(B))>tol
  error('This basis type is not fftable');
end

A = round(A);
B = round(B);

chebdata = jac.jfft.chebfft_overhead(N,opt);
if strcmpi(opt.normalization,'normal')
  C = jac.connection.integer_separation_connection_matrix(N,opt.alpha,opt.beta,A,B);
else
  error('Unrecognized polynomial normalization specification');
end

[data.A,data.B,data.N,data.C,data.chebdata] = deal(A,B,N,C,chebdata);
