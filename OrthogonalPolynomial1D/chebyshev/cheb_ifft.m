function[f] = cheb_ifft(F);

% [f] = cheb_ifft(F)
% Uses the IFFT to compute a modal to nodal expansion for the normalized
% Chebyshev polynomials
%
% Assumes the nodal evaluations f are sought at the Chebyhsev-Gauss nodal points.
% F must be a column vector (but this is vectorized across columns)
%
% 20080812 -- acn

N = size(F,1);
shift = (1+1/(2*N))*(0:(-1):-(N-1)).';
shift = 1/sqrt(2*pi)*exp(-i*pi*shift);
shift(1) = sqrt(2)*shift(1);

f = spdiags(shift,0,N,N)*F;
f = [0; real(flipud(f(2:end))); real(f)] + i*[0; -imag(flipud(f(2:end))); imag(f)];
f = ifft(ifftshift(f,1),[],1)*2*N;
f = real(f(1:N));
