function[data] = integer_separation_connection_overhead(N,gamma,delta,G,D)
% [data] = integer_separation_connection_overhead(N,gamma,delta,G,D)
%
%     Performs overhead computations for Fourier integer-parameter connections
%     and stores the information in data.

global handles;
jac = handles.speclab.orthopoly1d.jacobi;
fourier = handles.speclab.fourier;
sc_expand = fourier.connection.sc_expand;
sc_collapse = fourier.connection.sc_collapse;

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
C_even = jac.connection.integer_separation_connection_matrix(Nc,alpha,beta,A,B);
C_odd = jac.connection.integer_separation_connection_matrix(Nc-1,alpha+1,beta+1,A,B);

[data.C_even, data.C_odd, data.N, data.no_connect, data.GD] = ...
  deal(C_even, C_odd,N, no_connect, GD);

% modes = sc_expand(cmodes,smodes,N);
