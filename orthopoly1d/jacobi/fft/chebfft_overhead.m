function[data] = chebfft_overhead(N,varargin);
% [data] = chebfft_overhead(N,{normalization='normal', scale=1, points = 'gq'})
%
%     Computes overhead data for Chebyshev FFT usage.
%
%     Assumes the nodal evaluations f are located at the Chebyhsev-Gauss nodal
%     points. If they're located at different locations, make the optional input
%     points either 'grq' for Gauss-Radau or 'glq' for Gauss-Lobatto.

global handles;
opt = handles.common.input_schema({'normalization','scale', 'points'}, ...
  {'normal',1, 'gq'},[], varargin{:});

if strcmpi(opt.points, 'gq')
  %%%%%%%%%%% The FFT way -- almost always slower
  %shift = (1+1/(2*N))*(0:(-1):-(N-1)).';
  %shift = sqrt(2*pi)*exp(i*pi*shift);
  %shift(1) = shift(1)/sqrt(2);
  %shift = shift/(2*N);
  %ft = complex(zeros([2*N 1]));

  %[data.N, data.shift,data.points, data.ft] = ...
  %  deal(N,shift,opt.points,ft);
  %%%%%%%%%%%

  %%%%%%%%%%% The DCT way -- usually faster
  shift = (-1).^([0:(N-1)].')/sqrt(N/pi);;
  [data.N, data.shift,data.points] = deal(N,shift,opt.points);
  %%%%%%%%%%%

elseif strcmpi(opt.points, 'grq')
  error('Not yet implemented')
elseif strcmpi(opt.points, 'glq')
  error('Not yet implemented')
else
  error('Unrecognized flag for location of nodes')
end
