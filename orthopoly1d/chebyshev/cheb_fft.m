function[F] = cheb_fft(f);

% [F] = cheb_fft(f)
% Uses the FFT to compute a nodal to modal expansion for the normalized
% Chebyshev polynomials
%
% Assumes the nodal evaluations f are located at the Chebyhsev-Gauss nodal points.
% f must be a column vector (but this is vectorized across columns)
%
% 20080812 -- acn

N = size(f,1);

ft = [f; flipud(f)];

ft = fftshift(fft(ft,[],1),1)/(2*N);

ft = ft((N+1):(2*N),:);
shift = (1+1/(2*N))*(0:(-1):-(N-1)).';

shift = sqrt(2*pi)*exp(i*pi*shift);
shift(1) = shift(1)/sqrt(2);

ft = spdiags(shift,0,N,N)*ft;

F = real(ft);
