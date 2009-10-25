function[f] = jifft(F,varargin)
% [f] = jifft(F,{points='gq', alpha=-1/2, beta=-1/2, normalization='normal', scale=1})
%
%     If 2*alpha and 2*beta are odd, this function computes the nodal
%     evaluations from the Jacobi Polynomial modal coefficients F using the
%     IFFT. The locations of the nodes must be a **Chebyshev** fft-compatible
%     set of points. See chebifft. 

global packages;
jac = packages.speclab.orthopoly1d.jacobi;
la = packages.labtools.linalg;
inputs = {'points', 'alpha', 'beta', 'normalization', 'scale'};
defaults = {'gq', -1/2, -1/2, 'normal', 1};
opt = packages.labtools.input_schema(inputs, defaults, [], varargin{:});

[tf,A,B] = jac.fft.fftable(opt);
if not(tf)
  error('Cannot use the FFT when 2*alpha and 2*beta are not odd integers');
end
N = size(F,1);

C = jac.connection.integer_separation_connection_matrix(N,-1/2,-1/2,A,B);
F = la.triu_sparse_invert(C,F,'bandwidth',A+B+1);
f = jac.fft.chebifft(F,opt);
