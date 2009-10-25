function[data] = wfft_galerkin_overhead(N,varargin)
% [F] = wfft_galerkin_overhead(N,{s=1, t=0, scale=1, shift=0})
%
%     Computes overhead data required to perform the FFT 'galerkin' algorithm 
%     for the Wiener basis functions. See wfft_galerkin.m

global packages;
wiener = packages.speclab.wiener;
opt = wiener.defaults(varargin{:});
fourier = packages.speclab.fourier;

[fftable,S,T] = wiener.fft.fftable(opt);
if not(fftable)
  error('This basis set is not fft-able: s and t must be integers');
end

% Compute fft overhead
data.fftdata = fourier.fft.ffft_overhead(N);

% Compute connection overhead
data.connection = fourier.connection.integer_separation_connection_overhead(N,0,0,S,T);

[data.S, data.T, data.ST, data.opt] = deal(S,T,S+T,opt);
