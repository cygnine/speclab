function[modes] = sc_expand(cmodes,smodes,N)
% [modes] = sc_expand(cmodes,smodes,N)
%
%     Function helper: expands the even cosine-like modes cmodes and the odd
%     sine-like modes smodes into a length-N vector modes.

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
