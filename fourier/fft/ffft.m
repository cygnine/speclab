function[modes] = ffft(nodes,varargin);
% [modes] = ffft(nodes,{gamma=0,delta=0,shift=0,scale=1});
% 
%     The "Fourier FFT". ("fft" conflicts Matlab's builtin)
%
%     A wrapper for Matlab's FFT. This function scales the modes so that they
%     correspond to modal coefficients in the basis expansion
%     1/sqrt(2*pi)*exp(i*k*theta), where k is an integer index, negatively
%     biased when N = length(nodes) is even. If the inputs gamma or delta are
%     greater than 0 and are natural numbers, this function also performs the
%     connection problem so that the output is the modes for the (gamma,delta)
%     generalized Fourier expansion.

persistent input_schema connection integer_range
if isempty(connection)
  from labtools import input_schema
  from speclab.fourier.connection import positive_integer_separation_connection as connection
  from speclab.common import integer_range
end

opt = input_schema({'gamma','delta','shift','scale'},{0,0,0,1},[],varargin{:});
N = length(nodes);

% This commented-out way is the way I *should* be doing it....
% modes = fftshift(fft(fftshift(nodes)));

modes = fft(nodes);
modes = fftshift(modes)/N*sqrt(2*pi);

ks = integer_range(N);
phase = exp(-i*ks*pi/N);
phase(ks==0) = 1;
phase(mod(ks,2)==1) = phase(mod(ks,2)==1)*-1;
modes = phase.*modes;

if (opt.gamma+opt.delta)>0
  modes = connection(modes,0,0,opt.gamma,opt.delta);
end
