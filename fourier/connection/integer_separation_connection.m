function[modes] = integer_separation_connection(modes,gamma,delta,G,D)
% [MODES] = INTEGER_SEPARATION_CONNECTION(MODES,GAMMA,DELTA,G,D)
%
%     Connects the modal coefficients MODES from an unweighted Fourier expansion
%     of class (GAMMA,DELTA) to the modal coefficients for an expansion of class
%     (GAMMA+G,DELTA+D). The connection coefficients are not explicity computed. 

global handles;
pos_conn = handles.speclab.fourier.connection.positive_integer_separation_connection;
neg_conn = handles.speclab.fourier.connection.negative_integer_separation_connection;

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

modes = pos_conn(modes,gamma,delta,G_pos,D_pos);

G_neg = abs(min([0,G]));
D_neg = abs(min([0,D]));

modes = neg_conn(modes,gamma+G_pos,delta+D_pos,G_neg,D_neg);
