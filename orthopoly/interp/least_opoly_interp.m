function[fz] = least_opoly_interp(theta, f, z, varargin);
% least_opoly_interp -- Multidimensional interpolation via the (opoly) least interpolant
%
% fz = least_opoly_interp(theta, f, z, {basis=[], domain=[],  maxN=1e6/length(c)})
%
%     Using the 'orthogonal polynomial least interpolant' of de Boor and Ron,
%     this function evaluates the interpolant at the locations z, where the
%     interpolant is specified by the data (theta, f). The default is Legendre
%     polynomials on [x1, y1] x [x2 y2] x [x3 y3] x [xd, yd]. If the d x 2 array
%     'domain' is given, the values of xi, yi are extracted from this. If not,
%     then the domain is defined from the nodal positions:
%
%       x1 = min(theta(:,1))     y1 = max(theta(:,1))
%       x2 = min(theta(:,2))     y2 = max(theta(:,2))
%       x3 = min(theta(:,3))     y3 = max(theta(:,3))
%          .                        .
%          .                        .
%          .                        .
%
%     If however a Basis object is given as the optional input 'basis', then
%     that basis is assumed to be an appropriate orthogonal polynomial basis and
%     its evaluations are used, ignoring anything given in 'domain'.
% 
%     If z is very large, memory may become an issue for evaluating the
%     interpolant since a large Vandermonde matrix is then required. The
%     optional input maxN determines how many points to generate a Vandermonde
%     matrix for, and then a for-loop is employed to loop over all the points.

persistent least_opoly_lu least_opoly_coeffs myip_full myip space_dim subspace_dim strict_inputs
if isempty(least_opoly_lu)
  from labtools import strict_inputs
  from speclab.orthopoly.interp import least_opoly_lu least_opoly_coeffs
  from speclab.common.tensor import space_dimension as space_dim
  from speclab.common.tensor import subspace_dimension as subspace_dim
  myip_full = @(d,k) speye(space_dim(d,k));
  myip = @(d,k) speye(subspace_dim(d,k));
end

opt = strict_inputs({'basis', 'domain', 'maxN'}, {[], [], floor(1e6/size(theta,1))}, [], varargin{:});

if isempty(z)
  fz = [];
  return;
end

% First find domain for basis
if isempty(opt.basis)
  dim = size(theta,2);
  domain = zeros([dim 2]);
  bases = cell([dim 1]);
  for q = 1:dim
    if isempty(opt.domain)
      interval = [min(theta(:,q)) max(theta(:,q))];
    else
      interval = opt.domain(q,:);
    end
    bases{q} = LegendrePolynomialBasis('domain', interval);
  end
  TL = TensorProductBasis(bases{:});
  evalbasis = @(x,n) TL.evaluate(x,n+1);
else
  evalbasis = @(x,n) opt.basis.evaluate(x,n+1);
  dim = opt.basis.dim;
end

[l,u,p,v,k_storage] = least_opoly_lu(theta, 'basis', evalbasis, 'ip', myip_full);
c = least_opoly_coeffs(l,u,p,v,k_storage,dim,f, 'ip', myip);

maxN = opt.maxN;

Nz = size(z, 1);
fz = zeros([Nz 1]);
for q = 1:ceil(Nz/maxN)
  i1 = (q-1)*maxN + 1;
  i2 = q*maxN;
  if q == ceil(Nz/maxN)
    i2 = Nz;
  end

  V = evalbasis(z(i1:i2,:), 0:(length(c)-1));
  fz(i1:i2) = V*c;
end
