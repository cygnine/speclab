function[data] = wfft_collocation_overhead(N,varargin)
% [F] = wfft_collocation_overhead(f,{s=1, t=0, scale=1, shift=0})
%
%     Precomputes data necessary to implement the FFT collocation method for the
%     Wiener rational basis functions on N points. The nodal locations are
%     assumed to be the canonical Fourier points mapped to the real line. Both
%     parameters s and t must be integers.

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
[x,w] = wiener.quad.gauss_quadrature(N,wopt);
weight = wiener.weights.phase_shifted_sqrt_weight(x,opt);
fopt = struct('gamma', S, 'delta', T);
fourier_data = fourier.fft.ffft_overhead(N,fopt);

[data.weight, data.fourier_data] = deal(weight, fourier_data);
