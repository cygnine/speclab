function[modes] = negative_integer_separation_connection_online(modes,data)
% [modes] = negative_integer_separation_connection_online(modes,data)
%
%     Connects the modal coefficients modes from an unweighted Fourier expansion
%     of class (gamma,delta) to the modal coefficients for an expansion of class
%     (gamma-G,delta-D) for G,D>0. The values of G, D, gamma, and delta are
%     already set in the construction of the input data.

persistent sc_expand sc_collapse matinv
if isempty(sc_expand)
  from speclab.fourier.connection import sc_expand sc_collapse
  from labtools.linalg import triu_sparse_invert as matinv
end

if data.no_connect
  return
end

[cmodes,smodes] = sc_collapse(modes,data.N);

cmodes = matinv(data.C_even,cmodes,'bandwidth',data.GD+1);
smodes = matinv(data.C_odd,smodes,'bandwidth',data.GD+1);

modes = sc_expand(cmodes,smodes,data.N);
