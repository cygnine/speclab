function[F] = jfft(f,varargin)
% [F] = jfft(f,{points='gq', alpha=-1/2, beta=-1/2, normalization='normal', scale=1})
%
%     If 2*alpha and 2*beta are odd, this function computes the Jacobi
%     Polynomial modal coefficients of the function point-values f. The
%     locations of the nodes must be a **Chebyshev** fft-compatible set of
%     points. See chebfft. 

persistent input_schema connection chebfft fftable
if isempty(input_schema)
  from labtools import input_schema
  from speclab.orthopoly1d.jacobi.connection import integer_separation_connection_matrix as connection
  from speclab.orthopoly1d.jacobi.fft import chebfft as cfft
  from speclab.orthopoly1d.jacobi.fft import chebfft fftable
end

inputs = {'points', 'alpha', 'beta', 'normalization', 'scale'};
defaults = {'gq', -1/2, -1/2, 'normal', 1};
opt = input_schema(inputs, defaults, [], varargin{:});

tol = 1e-12;
[tf,A,B] = fftable(opt);
if not(tf)
  error('Cannot use the FFT when 2*alpha and 2*beta are not odd integers');
end
N = size(f,1);

C = connection(N,-1/2,-1/2,A,B);
F = chebfft(f,opt);
F = C*F;
