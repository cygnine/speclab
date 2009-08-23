%%% Example script for Wiener approximation
%%% wfft: Usage of the FFT for Wiener approximations.

% As with Jacobi/Fourier expansions, there are online/overhead and direct FFT
% methods. In addition, there are two algorithms for using the FFT: the
% `Galerkin' and the 'collocation' method. The Galerkin method is slower, but
% has less of a numerical penalty if you try to expand something decaying like
% 1/|x|^2 in s=12 functions.

% A very important caveat for the FFT for the Wiener expansions is that the
% terminal s modes (positive/negative) are garbage. They should be discarded.
% This is true for both the collocation and Galerkin methods. Over-integration
% by 2*s is more or less required.

% In general, the collocation and Galerkin methods are double-edged swords:
% - the collocation method is much better-suited for expansions where the modal
%   coefficients decay quickly.
% - the Galerkin method is much better-suited for expansions where the modal
%   coefficients decay slowly.
% If s is large (> ~10), the Galerkin method also suffers more due to roundoff
% error

clear
global handles;
wiener = handles.speclab.wiener;
wfft = wiener.fft;
irange = handles.speclab.common.integer_range;

opt.s = ceil(1 + 9*rand); % Just some natural number >= 2
% Now that s is determined, let's add on the over-integration penalty. This is
% not necessary for direct s=1 quadrature to get results that are not clearly
% wrong, but the s=1 quadrature rule requires 2*s over-integration to accurately
% compute the mass matrix anyway.
N = 400; Nq = N + 2*opt.s;

f = @(x) exp(-x.^2);
% For scaling: put 90% of the resolution inside the 1/e^4 ~ 0.018 window. 1/e^4
% for the Gaussian happens at x= +/- 2
% Unlike the finite-interval affine_scaling functions, N is a required input for
% infinite-interval scalings
temp = wiener.affine_scaling([-2,2], N, 's', opt.s, 'resolution_fraction', 0.9);
opt.scale = temp.scale; opt.shift= temp.shift;
% Need s=1 nodes for the fft: 
opt1 = opt; opt1.s = 1;

ks = irange(N);
[x,w] = wiener.quad.pi_gauss_quadrature(N,opt);
ws = wiener.eval.wiener_function(x,ks,opt);
direct_modes = ws'*(f(x).*w);


% This is a little subtle: for the FFT to work, we need the s=1 nodes. However,
% the scaling rubric computed earlier doesn't mean anything for the s=1 nodes.
% (It is *not* "90% of the nodes are inside [-2,2]", especially since we're
% using Nq here when we scaled for N.) Nevertheless, the same scaling parameter
% must be used for the FFT to give correct results.
[x,w] = wiener.quad.pi_gauss_quadrature(Nq,opt1);
% But must tell the fft which value of s we want:
fft_collocation_modes = wfft.wfft_collocation(f(x),opt);
fft_galerkin_modes = wfft.wfft_galerkin(f(x),opt);

% Online/overhead methods:
fftdata_collocation = wfft.wfft_collocation_overhead(Nq,opt);
fftdata_galerkin = wfft.wfft_galerkin_overhead(Nq,opt);

fft_collocation_online_modes = wfft.wfft_collocation_online(f(x),fftdata_collocation);
fft_galerkin_online_modes = wfft.wfft_galerkin_online(f(x),fftdata_galerkin);

s = opt.s;
% Discard garbage modes:
fft_collocation_modes = fft_collocation_modes(s+1:end-s);
fft_galerkin_modes = fft_galerkin_modes(s+1:end-s);
fft_collocation_online_modes = fft_collocation_online_modes(s+1:end-s);
fft_galerkin_online_modes = fft_galerkin_online_modes(s+1:end-s);

% You can now check the difference between the modes, e.g.:
err = norm(direct_modes - fft_collocation_modes);
fprintf('The difference between direct and collocation is %1.4e\n', err);
err = norm(direct_modes - fft_galerkin_modes);
fprintf('The difference between direct and Galerkin is %1.4e\n', err);

% And you can also check that the online/straightforward fft methods produce
% (nearly) identical results.

% As with the Fourier/Jacobi case, the routines
% - iwfft_collocation
% - iwfft_collocation_online
% - iwfft_galerkin
% - iwfft_galerkin_online
% exist and do what it looks like they do. However, there are no iwfft_overhead
% routines. Instead, use wfft_collocation_overhead and wfft_galerkin_overhead to
% generate the fftdata and feed that straight into the iwfft online routines. 
