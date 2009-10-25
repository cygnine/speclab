function[f] = iwfft_collocation(F,varargin)
% [f] = iwfft_collocation(F,{s=1, t=0, scale=1, shift=0})
%
%     Uses the IFFT 'collocation' approach to compute nodal evaluations for the
%     Wiener basis functions. The nodal locations are assumed to be the
%     canonical Fourier points mapped to the real line. Both parameters s and t
%     must be integers.

global packages;
wiener = packages.speclab.wiener;
fourier = packages.speclab.fourier;
opt = wiener.defaults(varargin{:});

[fftable,S,T] = wiener.fft.fftable(opt);
if not(fftable)
  error('This basis set is not fft-able: s and t must be integers');
end

wopt = opt;
wopt.s = 1; wopt.t = 0;
[x,w] = wiener.quad.gauss_quadrature(length(F),wopt);
weight = wiener.weights.phase_shifted_sqrt_weight(x,opt);
fopt = struct('gamma', S, 'delta', T);

f = fourier.fft.iffft(F,fopt);
f = f.*weight;
