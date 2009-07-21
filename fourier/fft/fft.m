function[modes] = fft(nodes,varargin);
% [MODES] = FFT(NODES,G,D);
%
%     A wrapper for Matlab's FFT. This function scales the modes so that they
%     correspond to modal coefficients in the basis expansion
%     1/sqrt(2*pi)*exp(i*k*theta), where k is an integer index, negatively
%     biased when N = length(modes) is even. If the inputs G or D are greater
%     than 0 and are natural numbers, this function also performs the connection
%     problem so that the output is the modes for the (G,D) generalized Fourier
%     expansion.

global handles;
opt = handles.common.InputSchema({'G','D'},{0,0},[],varargin{:});

fourier = handles.speclab.fourier;

modes = fft(nodes);
modes = fftshift(modes)/length(modes)*sqrt(2*pi);

if (opt.G+opt.D)>0
  modes = fourier.connection.integer_separation_connection(modes,0,0,opt.G,opt.D);
end
