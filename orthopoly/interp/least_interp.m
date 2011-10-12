function[fz] = least_interp(theta, f, z, varargin);
% least_interp -- Multidimensional interpolation via the least interpolant
%
% fz = least_interp(theta, f, z, [maxN=1e6/length(c)])
%
%     Using the 'least interpolant' of de Boor and Ron, this function evaluates
%     the interpolant at the locations z, where the interpolant is specified by
%     the data (theta, f).
%
%     If z is very large, memory may become an issue for evaluating the
%     interpolant since a large Vandermonde matrix is then required. The
%     optional input maxN determines how many points to generate a Vandermonde
%     matrix for, and then a for-loop is employed to loop over all the points.

persistent least_lu least_coeffs multimonomial
if isempty(least_lu)
  from speclab.orthopoly.interp import least_lu least_coeffs
  from speclab.monomials import multimonomial
end

if isempty(z)
  fz = [];
  return;
end

dim = size(theta,2);
[l,u,p,v,k_storage] = least_lu(theta);

c = least_coeffs(l,u,p,v,k_storage,dim,f);

if nargin>3
  maxN = varargin{1};
else
  maxN = floor(1e6/length(c));
end

Nz = size(z, 1);
fz = zeros([Nz 1]);
for q = 1:ceil(Nz/maxN)
  i1 = (q-1)*maxN + 1;
  i2 = q*maxN;
  if q == ceil(Nz/maxN)
    i2 = Nz;
  end

  V = multimonomial(z(i1:i2,:), 0:(length(c)-1), 'dim', dim);
  fz(i1:i2) = V*c;
end
