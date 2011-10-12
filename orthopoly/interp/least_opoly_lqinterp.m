function[fz] = least_opoly_lqinterp(theta, f, z, varargin);
% least_opoly_lqinterp -- Multidimensional interpolation via the (opoly) least interpolant
%
% fz = least_opoly_lqinterp(theta, f, z, {basis=[], domain=[],  maxN=1e6/length(c)})
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

persistent least_opoly_lu least_opoly_coeffs myip_full myip space_dim subspace_dim parser
if isempty(least_opoly_lu)
  from labtools import input_parser
  from speclab.orthopoly.interp import least_opoly_lqlu least_opoly_lqcoeffs
  from speclab.common.tensor import subspace_dimension as subspace_dim

  myip = @(d,k) speye(subspace_dim(d,k));

  [opt, parser] = input_parser({'basis', 'domain', 'maxN'}, ...
                               {[], [], floor(1e6/size(theta,1))}, ...
                               [], ...
                               varargin{:});
else
  parser.parse(varargin{:});
  opt = parser.Results;
end

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
      imean = mean(interval);
      interval = 1.2*(interval - imean) + imean;
    else
      interval = opt.domain(q,:);
    end
    bases{q} = LegendrePolynomialBasis('domain', interval);
  end
  TL = TensorProductBasis(bases, 'indexing', 1);
  evalbasis = @(x,n) TL.evaluate(x,n);
else
  evalbasis = @(x,n) opt.basis.evaluate(x,n);
  dim = opt.basis.dim;
end

[l,u,p,v,k_storage] = least_opoly_lqlu(theta, 'basis', evalbasis, 'ip', myip);
c = least_opoly_lqcoeffs(l,u,p,v,k_storage,dim,f, 'ip', myip);

maxN = opt.maxN;

Nz = size(z, 1);
fz = zeros([Nz 1]);
for q = 1:ceil(Nz/maxN)
  i1 = (q-1)*maxN + 1;
  i2 = q*maxN;
  if q == ceil(Nz/maxN)
    i2 = Nz;
  end

  V = evalbasis(z(i1:i2,:), 1:length(c));
  fz(i1:i2) = V*c;
end
