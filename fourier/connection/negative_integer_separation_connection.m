function[modes] = negative_integer_separation_connection(modes,gamma,delta,G,D)
% [modes] = negative_integer_separation_connection(modes,gamma,delta,G,D)
%
%     Connects the modal coefficients modes from an unweighted Fourier expansion
%     of class (gamma,delta) to the modal coefficients for an expansion of class
%     (gamma-G,delta-D) for G,D>0. The connection coefficients are not explicity
%     computed. 

persistent sc_collapse sc_expand connection matinv
if isempty(sc_collapse)
  from speclab.fourier.connection import sc_collapse sc_expand
  from speclab.orthopoly.jacobi.connection import integer_separation_connection_matrix as connection
  from labtools.linalg import triu_sparse_invert as matinv
end

if (G+D)<=0
  return
end

alpha = delta-1/2;
A = D;
beta = gamma-1/2;
B = G;
N = length(modes);

[cmodes,smodes] = sc_collapse(modes,N);
Nc = ceil((N+1)/2);
C_even = connection(Nc,alpha-A,beta-B,A,B);
C_odd = connection(Nc-1,alpha+1-A,beta+1-B,A,B);

cmodes = matinv(C_even,cmodes,'bandwidth',G+D+1);
smodes = matinv(C_odd,smodes,'bandwidth',G+D+1);

modes = sc_expand(cmodes,smodes,N);
