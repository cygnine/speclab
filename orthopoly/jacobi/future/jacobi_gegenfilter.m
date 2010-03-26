function[ft] = jacobi_gegenfilter(t,modes);

% [ft] = jacobi_gegenfilter(t,modes)
% Uses the Gibbs phenomenon removal procedure of Gottlieb/Shu to re-expand the
% Chebyshev expansion given by modes into a Gegenbauer reconstruction of a given
% order. This expansion is truncated and interpolated at the points t to produce
% the output ft. It is assumed that the points t are in the interval [-1,1].

global common;
prevpath = addpaths(common.bases.d1.fourier);

N = length(modes);
if mod(N,2)==0;
  Ns = (-(N/2):(N/2-1)).';
else
  Ns = -((N-1)/2):((N-1)/2);
end

m = min([round(N/10) 15]);

% For now use interpolation, but we know the exact modal expansion of e^(i*k*x),
% so we should really use that instead.
jpolyN = N+2*m;
modes = [modes; zeros([2*m 1])];
%fr = cheb_ifft(modes);
%[r,w] = jacobipoly_gq(jpolyN,-1/2,-1/2);

% Find nodal values by interpolating to r points
%fr = Psin(pi*r,Ns,0,0)*modes;

%jacobimodes = jacobi_rmatrix(jpolyN,-1/2,-1/2,m,m)*cheb_fft(fr);
jacobimodes = jacobi_rmatrix(jpolyN,-1/2,-1/2,m,m)*modes;
%jacobimodes = jacobipolyn(r,0:(jpolyN-1),-1/2,-1/2)'*spdiags(w,0,jpolyN,jpolyN)*fr;
jacobimodes((floor(m)+1):end) = 0;

ft = jacobipolyn(t,0:(jpolyN-1),m-1/2,m-1/2)*jacobimodes;

path(prevpath);
