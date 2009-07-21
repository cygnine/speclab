function[modes] = integer_separation_connection(modes,gamma,delta,G,D)
% [MODES] = INTEGER_SEPARATION_CONNECTION(MODES,GAMMA,DELTA,G,D)
%
%     Connects the modal coefficients MODES from an unweighted Fourier expansion
%     of class (GAMMA,DELTA) to the modal coefficients for an expansion of class
%     (GAMMA+G,DELTA+D). The connection coefficients are not explicity computed. 

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;

% Test if G,D are integers:
tol = 1e-8;
if (abs(round(G)-G)>tol) || (abs(round(D)-D)>tol)
  error('Inputs G and D must be integers');
end
G = round(G);
D = round(D);

% Force column vector
modes = modes(:);

if (G+D)==0
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

end

function[cmodes,smodes] = sc_collapse(modes,N)
% [CMODES,SMODES] = SC_COLLAPSE(MODES,N)
%
%     Subfunction helper: collapses the fourier MODES into even cosine-like
%     modes CMODES and odd sine-like modes SMODES.

N_even = mod(N,2)==0;

if N_even
  modes(1) = modes(1)/2;
  modes(end+1) = conj(modes(1));
end

% Matlab default implicit casting: double + boolean = double
N = N + N_even;

N_positive = floor(N/2);

cmodes = zeros([N_positive+1,1]);
smodes = zeros([N_positive,1]);

cmodes = flipud(modes(1:N_positive+1));
smodes = -cmodes(2:end);

cmodes(2:end) = cmodes(2:end) + modes(N_positive+2:end);
smodes = smodes + modes(N_positive+2:end);
cmodes(1) = cmodes(1)*sqrt(2);

cmodes = cmodes/2;
smodes = smodes/2;

end

function[modes] = sc_expand(cmodes,smodes,N)
% [MODES] = SC_EXPAND(CMODES,SMODES,N)
%
%     Subfunction helper: expands the even cosine-like modes CMODES and the odd
%     sine-like modes SMODES into a length-N vector MODES.

N_even = mod(N,2)==0;
n = floor(N/2);

pmodes = cmodes;
nmodes = smodes;

% 'positive' modes = cmodes
pmodes(1) = cmodes(1)*sqrt(2);
pmodes(2:end) = cmodes(2:end) + smodes;

% 'negative' modes = smodes
nmodes = cmodes(2:end) - smodes;
nmodes = flipud(nmodes);

% Splice them together
modes = [nmodes;pmodes];
if N_even
  modes(1) = modes(1) + conj(modes(end));
  modes = modes(1:end-1);
end

end
