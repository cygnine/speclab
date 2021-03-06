function[data] = integer_separation_connection_overhead(N,gamma,delta,G,D)
% [data] = integer_separation_connection_overhead(N,gamma,delta,G,D)
%
%     Performs overhead computations for Fourier integer-parameter connections
%     and stores the information in data.

persistent connection
if isempty(connection)
  from speclab.orthopoly.jacobi.connection import integer_separation_connection_matrix as connection
end

GD = G+D;
no_connect = GD==0;

alpha = delta-1/2;
A = D;
beta = gamma-1/2;
B = G;

% At a future time: should get online/overhead routines for
% sc_expand/collapse.
%[cmodes,smodes] = sc_collapse(modes,N);

Nc = ceil((N+1)/2);
C_even = connection(Nc,alpha,beta,A,B);
C_odd = connection(Nc-1,alpha+1,beta+1,A,B);

[data.C_even, data.C_odd, data.N, data.no_connect, data.GD] = ...
  deal(C_even, C_odd,N, no_connect, GD);

% modes = sc_expand(cmodes,smodes,N);
