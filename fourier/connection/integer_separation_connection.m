function[modes] = integer_separation_connection(modes,gamma,delta,G,D)
% [modes] = integer_separation_connection(modes,gamma,delta,G,D)
%
%     Connects the modal coefficients modes from an unweighted Fourier expansion
%     of class (gamma,delta) to the modal coefficients for an expansion of class
%     (gamma+G,delta+D). The connection coefficients are not explicity computed. 

persistent pconnection nconnection
if isempty(pconnection)
  from speclab.fourier.connection import positive_integer_separation_connection as pconnection
  from speclab.fourier.connection import nositive_integer_separation_connection as nconnection
end

% Test if G,D are integers:
tol = 1e-8;
if (abs(round(G)-G)>tol) || (abs(round(D)-D)>tol)
  error('Inputs G and D must be integers');
end
G = round(G);
D = round(D);

% Force column vector
modes = modes(:);

if (G==0) & (D==0)
  return 
end

G_pos = max([0,G]);
D_pos = max([0,D]);

modes = pconnection(modes,gamma,delta,G_pos,D_pos);

G_neg = abs(min([0,G]));
D_neg = abs(min([0,D]));

modes = nconnection(modes,gamma+G_pos,delta+D_pos,G_neg,D_neg);
