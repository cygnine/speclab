function[f] = iwfft_collocation(F,varargin)
% [f] = iwfft_collocation(F,{s=1, t=0, scale=1, shift=0})
%
%     Uses the IFFT 'collocation' approach to compute nodal evaluations for the
%     Wiener basis functions. The nodal locations are assumed to be the
%     canonical Fourier points mapped to the real line. Both parameters s and t
%     must be integers.

persistent defaults fftable gauss_quadrature wf iffft
if isempty(defaults)
  from speclab.wiener import defaults
  from speclab.wiener.fft import fftable
  from speclab.wiener.quad import gauss_quadrature
  from speclab.wiener.weights import phase_shifted_sqrt_weight as wf
  from speclab.fourier.fft import iffft
end

opt = defaults(varargin{:});

[tf,S,T] = fftable(opt);
if not(tf)
  error('This basis set is not fft-able: s and t must be integers');
end

wopt = opt;
wopt.s = 1; wopt.t = 0;
[x,w] = gauss_quadrature(length(F),wopt);
weight = wf(x,opt);
fopt = struct('gamma', S, 'delta', T);

f = iffft(F,fopt);
f = f.*weight;
