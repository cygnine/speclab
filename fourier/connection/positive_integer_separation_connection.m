function[modes] = positive_integer_separation_connection(modes,gamma,delta,G,D)
% [modes] = positive_integer_separation_connection(modes,gamma,delta,G,D)
%
%     Connects the modal coefficients modes from an unweighted Fourier expansion
%     of class (gamma,delta) to the modal coefficients for an expansion of class
%     (gamma+G,delta+D) for G,D>0. The connection coefficients are not explicity
%     computed. It is assumed that the inputs G and D are integers.

persisten jac fourier sc_expand sc_collapse
if isempty(jac)
  imp speclab.orthopoly1d.jacobi as jac
  from speclab import fourier
  from fourier.connection import sc_expand sc_collapse
end

%global packages;
%jac = packages.speclab.orthopoly1d.jacobi;
%fourier = packages.speclab.fourier;
%sc_expand = fourier.connection.sc_expand;
%sc_collapse = fourier.connection.sc_collapse;

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
