function[modes] = ffft(nodes,varargin);
% [MODES] = FFFT(NODES,{GAMMA=0,DELTA=0,SHIFT=0,SCALE=1});
%
%     A wrapper for Matlab's FFT. This function scales the modes so that they
%     correspond to modal coefficients in the basis expansion
%     1/sqrt(2*pi)*exp(i*k*theta), where k is an integer index, negatively
%     biased when N = length(modes) is even. If the inputs GAMMA or DELTA are
%     greater than 0 and are natural numbers, this function also performs the
%     connection problem so that the output is the modes for the (GAMMA,DELTA)
%     generalized Fourier expansion.

global handles;
opt = handles.common.InputSchema({'gamma','delta','shift','scale'},{0,0,0,1},[],varargin{:});
conn = handles.speclab.fourier.connection.positive_integer_separation_connection;
N = length(nodes);

%modes = fftshift(fft(fftshift(nodes)));
%phase = 

modes = fft(nodes);
modes = fftshift(modes)/N*sqrt(2*pi);

ks = handles.speclab.common.integer_range(N);
phase = exp(-i*ks*pi/N);
phase(ks==0) = 1;
phase(mod(ks,2)==1) = phase(mod(ks,2)==1)*-1;
modes = phase.*modes;

if (opt.gamma+opt.delta)>0
  modes = conn(modes,0,0,opt.gamma,opt.delta);
end
