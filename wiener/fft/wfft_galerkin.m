function[F] = wfft_galerkin(f,varargin)
% [F] = wfft_galerkin(f,{s=1, t=0, scale=1, shift=0})
%
%     Uses the FFT 'galerkin' approach to compute modal coefficients for the
%     Wiener basis functions. The nodal locations are assumed to be the
%     canonical Fourier points mapped to the real line. Both parameters s and t
%     must be integers.

global handles;
wiener = handles.speclab.wiener;
fourier = handles.speclab.fourier;
opt = wiener.defaults(varargin{:});
wconnect = fourier.connection.positive_integer_separation_connection;

[fftable,S,T] = wiener.fft.fftable(opt);
if not(fftable)
  error('This basis set is not fft-able: s and t must be integers');
end

% First compute (s=1,t=0) *unweighted* modes
F = fourier.fft.ffft(f);

% Then compute same modes for f/w^s
F = wiener.operators.wiener_weight_divide(F,opt);

% Finally connect these from s=1 ---> s=s modes
if S>0
  F = wconnect(F,0,0,S,T);
  F(end-S-T:end) = 0;
end

F(1:(S+T+1)) = 0;
