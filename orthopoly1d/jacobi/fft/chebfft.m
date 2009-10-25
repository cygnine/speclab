function[F] = chebfft(f,varargin);
% [F] = chebfft(f,{normalization='normal', scale=1, points = 'gq'})
%
%     Uses the FFT to compute a nodal to modal expansion for the normalized
%     Chebyshev polynomials
%
%     Assumes the nodal evaluations f are located at the Chebyhsev-Gauss nodal
%     points. If they're located at different locations, make the optional input
%     points either 'grq' for Gauss-Radau or 'glq' for Gauss-Lobatto.
%     f must be a column vector (but this is vectorized across columns)

global packages;
opt = packages.labtools.input_schema({'normalization','scale', 'points'}, ...
  {'normal',1, 'gq'},[], varargin{:});

if strcmpi(opt.points, 'gq')
  N = size(f,1);
  ft = [f; flipud(f)];

  ft = fftshift(fft(ft,[],1),1)/(2*N);

  ft = ft((N+1):(2*N),:);
  shift = (1+1/(2*N))*(0:(-1):-(N-1)).';

  shift = sqrt(2*pi)*exp(i*pi*shift);
  shift(1) = shift(1)/sqrt(2);

  ft = spdiags(shift,0,N,N)*ft;

  F = real(ft);
elseif strcmpi(opt.points, 'grq')
  error('Not yet implemented')
elseif strcmpi(opt.points, 'glq')
  error('Not yet implemented')
else
  error('Unrecognized flag for location of nodes')
end
