function[fz] = least_opoly_lqinterp(theta, f, z, varargin);
% least_opoly_lqinterp -- Multidimensional interpolation via the least
%                         orthogonal interpolant
%
% fz = least_opoly_lqinterp(theta, f, z, {basis='legendre', ip=eye(), maxN=1e6/length(c)})
%
%     Using the least orthogonal polynomial, this function evaluates the
%     interpolant at the locations z, where the interpolant is specified by the
%     data (theta, f). The default is tensor-product Legendre polynomials
%     defined using univariate intervals that tightly circumscribe the nodes.
%
%     OPTIONAL INPUTS:
%     basis is as in least_opoly_lqlu: either a TensorProductBasis instance, or
%     a function handle supporting appropriate syntax.
%
%     ip is a function handle returning Gram matrices for certain degrees: see
%     least_opoly_lqlu.
% 
%     maxN is a positive integer for doing batch computations if z is very
%     large: memory may become an issue for evaluating the interpolant since a
%     large Vandermonde matrix is then required. maxN determines how many
%     points to generate a Vandermonde matrix for, and then a for-loop is
%     employed to loop over all the points.

persistent least_opoly_lu myip space_dim subspace_dim parser
if isempty(least_opoly_lu)
  from labtools import input_parser
  from speclab.orthopoly.interp import least_opoly_lqlu 
  from speclab.common.tensor import polynomial_subspace_dimension as subspace_dim

  ip = @(d,k) speye(subspace_dim(d,k));

  [opt, parser] = input_parser({'basis', 'ip', 'maxN'}, ...
                               {[], [], floor(1e6/size(theta,1))}, ...
                               [], ...
                               varargin{:});
else
  parser.parse(varargin{:});
  opt = parser.Results;
end

% Trivial?
if isempty(z)
  fz = [];
  return;
end

if isempty(opt.ip)
  opt.ip = ip;
end

%% First find domain for basis
%if isempty(opt.basis)
%  dim = size(theta,2);
%  domain = zeros([dim 2]);
%  bases = cell([dim 1]);
%  for q = 1:dim
%    if isempty(opt.domain)
%      interval = [min(theta(:,q)) max(theta(:,q))];
%      imean = mean(interval);
%      interval = 1.2*(interval - imean) + imean;
%    else
%      interval = opt.domain(q,:);
%    end
%    bases{q} = LegendrePolynomialBasis('domain', interval);
%  end
%  TL = TensorProductBasis(bases, 'indexing', 1);
%  evalbasis = @(x,n) TL.evaluate(x,n);
%else
%  evalbasis = @(x,n) opt.basis.evaluate(x,n);
%  dim = opt.basis.dim;
%end

%[l,u,p,v,k_storage,H] = least_opoly_lqlu(theta, 'basis', evalbasis, 'ip', myip);
[l,u,p,H,v,k_storage] = least_opoly_lqlu(theta, 'basis', opt.basis, 'ip', opt.ip);
if isa(opt.basis, 'TensorProductBasis')
  f = f.*sqrt(opt.basis.weight(theta));
end
c = H'*inv(l*u)*p*f;

maxN = opt.maxN;

Nz = size(z, 1);
fz = zeros([Nz size(f,2)]);
for q = 1:ceil(Nz/maxN)
  i1 = (q-1)*maxN + 1;
  i2 = q*maxN;
  if q == ceil(Nz/maxN)
    i2 = Nz;
  end

  %V = evalbasis(z(i1:i2,:), 1:length(c));
  if isa(opt.basis, 'TensorProductBasis')
    V = opt.basis(z(i1:i2,:), opt.basis.range(length(c)));
  else
    % Assume 1-based indexing
    V = opt.basis(z(i1:i2,:), 1:length(c));
  end
  fz(i1:i2,:) = V*c;
end
