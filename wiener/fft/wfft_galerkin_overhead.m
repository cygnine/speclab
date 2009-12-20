function[data] = wfft_galerkin_overhead(N,varargin)
% [F] = wfft_galerkin_overhead(N,{s=1, t=0, scale=1, shift=0})
%
%     Computes overhead data required to perform the FFT 'galerkin' algorithm 
%     for the Wiener basis functions. See wfft_galerkin.m

persistent defaults fftable ffft_overhead connection
if isempty(defaults)
  from speclab.wiener import defaults
  from speclab.wiener.fft import fftable
  from speclab.fourier.fft import ffft_overhead
  from speclab.fourier.connection import integer_separation_connection_overhead as connection
end

opt = defaults(varargin{:});

[tf,S,T] = fftable(opt);
if not(tf)
  error('This basis set is not fft-able: s and t must be integers');
end

% Compute fft overhead
data.fftdata = ffft_overhead(N);

% Compute connection overhead
data.connection = connection(N,0,0,S,T);

[data.S, data.T, data.ST, data.opt] = deal(S,T,S+T,opt);
