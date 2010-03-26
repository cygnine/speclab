function[f] = chebifft(F,varargin);
% [f] = chebifft(F,{normalization='normal', scale=1, points = 'gq'})
%
%     Uses the IFFT to compute the nodal evaluations f from the modal coefficients
%     F.
%
%     Assumes the nodal evaluations f are located at the Chebyhsev-Gauss nodal
%     points. If they're located at different locations, make the optional input
%     points either 'grq' for Gauss-Radau or 'glq' for Gauss-Lobatto.
%     F must be a column vector (but this is vectorized across columns)

persistent input_schema
if isempty(input_schema)
  from labtools import input_schema
end

opt = input_schema({'normalization','scale', 'points'}, ...
  {'normal',1, 'gq'},[], varargin{:});

if strcmpi(opt.points, 'gq')

  N = size(F,1);
  shift = (1+1/(2*N))*(0:(-1):-(N-1)).';
  shift = 1/sqrt(2*pi)*exp(-i*pi*shift);
  shift(1) = sqrt(2)*shift(1);

  f = spdiags(shift,0,N,N)*F;
  f = [0; real(flipud(f(2:end))); real(f)] + i*[0; -imag(flipud(f(2:end))); imag(f)];
  f = ifft(ifftshift(f,1),[],1)*2*N;
  f = real(f(1:N));

elseif strcmpi(opt.points, 'grq')
  error('Not yet implemented')
elseif strcmpi(opt.points, 'glq')
  error('Not yet implemented')
else
  error('Unrecognized flag for location of nodes')
end
