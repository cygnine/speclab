function[modes] = positive_integer_separation_connection(modes,gamma,delta,G,D)
% [MODES] = POSITIVE_INTEGER_SEPARATION_CONNECTION(MODES,GAMMA,DELTA,G,D)
%
%     Connects the modal coefficients MODES from an unweighted Fourier expansion
%     of class (GAMMA,DELTA) to the modal coefficients for an expansion of class
%     (GAMMA+G,DELTA+D) for G,D>0. The connection coefficients are not explicity
%     computed. 

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
fourier = handles.speclab.fourier;
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
C_even = jac.connection.integer_separation_connection_matrix(Nc,alpha,beta,A,B);
C_odd = jac.connection.integer_separation_connection_matrix(Nc-1,alpha+1,beta+1,A,B);

cmodes = C_even*cmodes;
smodes = C_odd*smodes;

modes = sc_expand(cmodes,smodes,N);
