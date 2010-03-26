function[C] = integer_separation_connection_matrix(N,alpha,A,varargin)
% [C] = integer_separation_connection_matrix(N,alpha,A,{normalization='normal'})
%
%     Returns the sparse N x N connection matrix C that connects modal
%     coefficients from Laguerre polynomials of class (alpha) to Laguerre
%     polynomials of class (alpha+A) when A is an integer. This connection is
%     sparse. Therefore, the matrix C is a sparse-matrix datatype. 
%
%     The optional input normalization defines the normalization of the
%     polynomials P, which affects the matrix C.

persistent input_schema x_times_p
if isempty(input_schema)
  from labtools import input_schema
  from speclab.orthopoly.laguerre import x_times_p
end
opt = input_schema({'normalization'}, {'normal'}, [],varargin{:});

C = spalloc(N,N,A*N - A*(A+1)/2);
C = C + speye(N);

ns = (0:N).';

for q = 1:A
  mu = x_times_p(ns,alpha+q,opt);
  C_temp = spdiags(mu,[0,-1],N,N).';
  C = C_temp*C;
end
