function[data] = jifft_overhead(N,varargin)
% [data] = jifft_overhead(N,{points='gq', alpha=-1/2, beta=-1/2, normalization='normal', scale=1})
%
%     If 2*alpha and 2*beta are odd, this function computes the overhead data
%     necessary for using the N-point Jacobi IFFT.

global handles;
jac = handles.speclab.orthopoly1d.jacobi;
inputs = {'points', 'alpha', 'beta', 'normalization', 'scale'};
defaults = {'gq', -1/2, -1/2, 'normal', 1};
opt = handles.common.InputSchema(inputs, defaults, [], varargin{:});

[tf,A,B] = jac.jfft.fftable(opt);

chebdata = jac.jfft.chebifft_overhead(N,opt);
if strcmpi(opt.normalization,'normal')
  C = jac.connection.integer_separation_connection_matrix(N,-1/2,-1/2,A,B);
else
  error('Unrecognized polynomial normalization specification');
end

[data.A,data.B,data.N,data.C,data.chebdata] = deal(A,B,N,C,chebdata);
