function[data] = wfft_collocation_overhead(N,varargin)
% [F] = wfft_collocation_overhead(f,{s=1, t=0, scale=1, shift=0})
%
%     Precomputes data necessary to implement the FFT collocation method for the
%     Wiener rational basis functions on N points. The nodal locations are
%     assumed to be the canonical Fourier points mapped to the real line. Both
%     parameters s and t must be integers.

persistent defaults fftable gauss_quadrature wf ffft_overhead
if isempty(defaults)
  from speclab.wiener import defaults
  from speclab.wiener.fft import fftable
  from speclab.wiener.quad import gauss_quadrature
  from speclab.wiener.weights import phase_shifted_sqrt_weight as wf
  from speclab.fourier.fft import ffft_overhead
end

opt = defaults(varargin{:});

[tf,S,T] = fftable(opt);
if not(tf)
  error('This basis set is not fft-able: s and t must be integers');
end

wopt = opt;
wopt.s = 1; wopt.t = 0;
[x,w] = gauss_quadrature(N,wopt);
weight = wf(x,opt);
fopt = struct('gamma', S, 'delta', T);
fourier_data = ffft_overhead(N,fopt);

[data.weight, data.fourier_data] = deal(weight, fourier_data);
