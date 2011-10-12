function[N] = marginalpoly_subspace_dim(dim, p)
% marginalpoly_subspace_dim -- dimension of subspace spanned by polynomials
%
% N = marginalpoly_subspace_dim(dim, p)
%
%     Returns the dimension of the polynomial space of dim-variate polynomials
%     of marginal degree exactly equal to p. Is vectorized in p.

psize = size(p);
p = p(:);
[p,ordering] = sort(p);

N = zeros(size(p));

if dim == 0
  N = zeros(size(p));
elseif dim == 1
  N = ones(size(p));
else
  P = max(p);
  Pinds = (0:P).';
  Ps = zeros(size(Pinds));
  for q = 0:dim-1
    Ps = Ps + nchoosek(dim, q)*(Pinds.^q);
  end

  pflags = [1; 1+find(diff(p)>0)];
  pflags = [pflags; length(p)];

  % Now distribute
  pmin = p(1);
  for q = 1:(length(pflags)-1)
    plevel = p(pflags(q));
    if plevel > -1
      N(pflags(q):(pflags(q+1))) = Ps(plevel+1);
    end
  end
end

% And reorder
N(ordering) = N;
N = reshape(N, psize);
