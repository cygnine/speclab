%%% Example script for Fourier approximation
%%% Introduction: calling syntax of basic routines

% Most of this is similar to the orthogonal polynomial case. Two major
% differences exist:
% - Modal coefficients indexing is not 0,1,2,3,..., but instead -N/2, -N/2+1,
% ..., N/2-1 if N is even, with a similar statement if N is odd. This causes
% complications in doing truncations. M-files do not yet exist for easy
% truncation of integer-indexed modes (e.g. Fourier and Wiener).
% - Derivatives are not implemented beyond first order.
%
% The standard interval of approximation for the Fourier series is [-pi, pi].

global handles;
fourier = handles.speclab.fourier;
irange = handles.speclab.common.integer_range;

% Quadrature+scaling doesn't really change
N = 100;
fopt = fourier.affine_scaling([-3,-1], 'resolution_fraction', 0.64, ...
                                       'N', N); % No reason, just for kicks
[theta,w] = fourier.quad.gauss_quadrature(N,fopt);

ks = irange(N);
% ks is a vector of modal coefficient indices. The goal is that if you tell me
% you want 7 modes, I'll return [-3, -2, -1, 0, 1, 2, 3]'. If you tell me you
% want 6 modes, I'll return [-3, -2, -1, 0, 1, 2]'. This "bias" towards
% negatives modes for N even is consistent throughout Speclab. The function
% handles.speclab.common.integer_range (or irange here) does this N <---> ks
% translation for you.
fs = fourier.eval.fseries(theta, ks, fopt);

fprintf('Using an N-point Gauss quadrature to compute NxN mass matrices when N is even...\n\n');
% The Fourier-gauss quadrature rules are degenerate for N even, meaning that
% they don't integrate 2*N functions exactly; rather only 2*N-1. Since for N
% even, there is a bias towards the negative mode (i.e. mode -N/2 is computed,
% but not +N/2), the (1,1) entry, effectively requiring exact integration of the
% 2*N'th function, will be incorrect. However, this is not the case for
% gamma=delta=0, which is the canonical Fourier case. Why not? Because this case
% is indeed special, and the quadrature rule is more accurate.
mass = fs'*spdiags(w,0,N,N)*fs;
fprintf('With gamma=delta=0:\n');
fprintf('The mass matrix differs from the identity by %1.3e\n\n', norm(mass-eye(N)));

% However, if we take, gamma, delta not 0, we'll see the (1,1) entry will not be
% correct.
fopt.gamma = -1/2 + 10*rand;
fopt.delta = -1/2 + 10*rand;
% Note that we're keeping the scaling parameters from gamma=delta=0, which are
% meaningless in this context other than they define some random interval of
% approximation..
[theta,w] = fourier.quad.gauss_quadrature(N,fopt);
fs = fourier.eval.fseries(theta,ks,fopt);
mass = fs'*spdiags(w,0,N,N)*fs;
str = ['With gamma=%1.3f and delta=%1.3f:\nThe mass matrix (2:N,2:N) differs '...
       'from the identity by %1.3e\n'];
fprintf(str, fopt.gamma, fopt.delta, norm(mass(2:N,2:N)-eye(N-1)));
fprintf('However, the (1,1) entry, mathematically 1, is computed as %1.3f\n', mass(1,1));

% If N is odd, there is no such degeneracy and the Gauss quadrature rule can
% exactly compute the mass matrix for any gamma,delta with a problem.
