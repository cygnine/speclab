function[data] = jfft_overhead(N,varargin)
% [data] = jfft_overhead(N,{points='gq', alpha=-1/2, beta=-1/2, normalization='normal', scale=1})
%
%     If 2*alpha and 2*beta are odd, this function computes the overhead data
%     necessary for using the N-point Jacobi FFT.

persistent input_schema fftable chebfft_overhead connection
if isempty(input_schema)
  from labtools import input_schema
  from speclab.orthopoly.jacobi.fft import fftable chebfft_overhead
  from speclab.orthopoly.jacobi.connection import integer_separation_connection_matrix as connection
end

inputs = {'points', 'alpha', 'beta', 'normalization', 'scale'};
defaults = {'gq', -1/2, -1/2, 'normal', 1};
opt = input_schema(inputs, defaults, [], varargin{:});

[tf,A,B] = fftable(opt);
if not(tf)
  error('Cannot use the FFT when 2*alpha and 2*beta are not odd integers');
end

chebdata = chebfft_overhead(N,opt);
if strcmpi(opt.normalization,'normal')
  C = connection(N,-1/2,-1/2,A,B);
else
  error('Unrecognized polynomial normalization specification');
end

[data.A,data.B,data.N,data.C,data.chebdata] = deal(A,B,N,C,chebdata);
