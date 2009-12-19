function[modes] = positive_integer_separation_connection_online(modes,data);
% [modes] = positive_integer_separation_connection_online(modes,data);
%
%     Connects the modal coefficients modes from an unweighted Fourier expansion
%     of class (gamma,delta) to the modal coefficients for an expansion of class
%     (gamma+G,delta+D) for G,D>0. The values of G, D, gamma, and delta were
%     already set in the construction of the input data. 

persistent sc_expand sc_collapse
if isempty(sc_expand)
  from speclab.fourier.connection import sc_expand sc_collapse
end

if data.no_connect
  return
end

[cmodes,smodes] = sc_collapse(modes,data.N);

cmodes = data.C_even*cmodes;
smodes = data.C_odd*smodes;

modes = sc_expand(cmodes,smodes,data.N);
