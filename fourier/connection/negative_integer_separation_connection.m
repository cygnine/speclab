function[modes] = negative_integer_separation_connection(modes,gamma,delta,G,D)
% [modes] = negative_integer_separation_connection(modes,gamma,delta,G,D)
%
%     Connects the modal coefficients modes from an unweighted Fourier expansion
%     of class (gamma,delta) to the modal coefficients for an expansion of class
%     (gamma-G,delta-D) for G,D>0. The connection coefficients are not explicity
%     computed. 

global packages;
jac = packages.speclab.orthopoly1d.jacobi;
fourier = packages.speclab.fourier;
la = packages.common.linalg;
sc_expand = fourier.connection.sc_expand;
sc_collapse = fourier.connection.sc_collapse;

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
C_even = jac.connection.integer_separation_connection_matrix(Nc,alpha-A,beta-B,A,B);
C_odd = jac.connection.integer_separation_connection_matrix(Nc-1,alpha+1-A,beta+1-B,A,B);

cmodes = la.triu_sparse_invert(C_even,cmodes,'bandwidth',G+D+1);
smodes = la.triu_sparse_invert(C_odd,smodes,'bandwidth',G+D+1);

modes = sc_expand(cmodes,smodes,N);
