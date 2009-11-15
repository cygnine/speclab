function[f] = iwfft_galerkin(F,varargin)
% [f] = iwfft_galerkin(F,{s=1, t=0, scale=1, shift=0})
%
%     Uses the IFFT 'galerkin' approach to compute function evaluations f from the
%     Wiener basis modal coefficients F. The nodal locations are assumed to be
%     the canonical Fourier points mapped to the real line. Both parameters s
%     and t must be integers.

persistent wiener fourier wconnect
if isempty(wiener)
  from speclab import wiener
  from speclab import fourier
  from speclab.fourier.connection import negative_integer_separation_connection as wconnect
end

%global packages;
%wiener = packages.speclab.wiener;
%fourier = packages.speclab.fourier;
%wconnect = fourier.connection.negative_integer_separation_connection;

opt = wiener.defaults(varargin{:});

[fftable,S,T] = wiener.fft.fftable(opt);
if not(fftable)
  error('This basis set is not fft-able: s and t must be integers');
end

% First connect down from unweighted s=s ----> s=1 modes
if S>0
  F = wconnect(F,S,T,S,T);
end

% Now compute (s=1,t=0) modes for f*w^s
F = wiener.operators.wiener_weight_multiply(F,opt);

% Finally compute nodal evaluations
f = fourier.fft.iffft(F);
