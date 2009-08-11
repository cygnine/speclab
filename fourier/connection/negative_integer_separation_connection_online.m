function[modes] = negative_integer_separation_connection_online(modes,data)
% [modes] = negative_integer_separation_connection_online(modes,data)
%
%     Connects the modal coefficients modes from an unweighted Fourier expansion
%     of class (gamma,delta) to the modal coefficients for an expansion of class
%     (gamma-G,delta-D) for G,D>0. The values of G, D, gamma, and delta are
%     already set in the construction of the input data.

if data.no_connect
  return
end

global handles;
fourier = handles.speclab.fourier;
la = handles.common.linalg;
sc_expand = fourier.connection.sc_expand;
sc_collapse = fourier.connection.sc_collapse;

[cmodes,smodes] = sc_collapse(modes,data.N);

cmodes = la.triu_sparse_invert(data.C_even,cmodes,'bandwidth',data.GD+1);
smodes = la.triu_sparse_invert(data.C_odd,smodes,'bandwidth',data.GD+1);

modes = sc_expand(cmodes,smodes,data.N);
