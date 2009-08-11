function[cmodes,smodes] = sc_collapse(modes,N)
% [cmodes,smodes] = sc_collapse(modes,N)
%
%     Function helper: collapses the fourier modes into even cosine-like
%     modes cmodes and odd sine-like modes smodes.

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
