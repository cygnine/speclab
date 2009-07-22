function[modes] = sc_expand(cmodes,smodes,N)
% [MODES] = SC_EXPAND(CMODES,SMODES,N)
%
%     Function helper: expands the even cosine-like modes CMODES and the odd
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
