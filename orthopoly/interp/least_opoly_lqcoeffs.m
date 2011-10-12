function[c] = least_opoly_lqcoeffs(l,u,p,v,k,D,f, varargin)
% least_opoly_lqcoeffs -- Computes basis coefficients for interpolation
%
% c = least_opoly_lqcoeffs(l,u,p,v,k,dim,f, [ip=alpha!, connection=speye])
%
%     Using the 5-tuple (l,u,p,v,k) that is output from de Boor's algorithm
%     deboor_lu, this function uses the point-evaluations f to compute the
%     basis coefficients necessary for interpolation. Obviously, f must
%     be an N-vector, where N is the number of rows of W (or l, or p, or u).
%     Unfortunately, the dimension of the problem cannot be divined from the
%     5-tuple, so the extra input dim denoting the dimension is necessary.
%
%     The interpolant may be evaluated as 
%
%       p(x) = sum_n c(n) * p(x,n),
%
%     where p(x,n) is a function that evaluates the n'th basis function at the
%     location x. (e.g. speclab.monomials.multimonomial is one such function). The
%     output c has length size(W,2). 
%
%     Note that this function is vectorized in the columns of f, meaning that if
%     f is a collection of column vectors, this function returns a collection of
%     column vectors containing the basis representation.

persistent invu invl dim subdim ip spdiag parser
if isempty(invu)
  from labtools import spdiag input_parser
  from labtools.linalg import triu_back_substitute as invu
  from labtools.linalg import tril_forward_substitute as invl

  from speclab.common.tensor import space_dimension as dim
  from speclab.common.tensor import subspace_dimension as subdim

  ip = @(d,k) speye(subdim(d,k));

  [opt, parser] = input_parser({'ip', 'connection'}, ...
                               {ip, @speye}, ...
                               [], ...
                               varargin{:});
else
  parser.parse(varargin{:});
  opt = parser.Results;
end

% These are the coefficients for the cardinal interpolants under the de Boor
% inner product. The remainder of the code translates these coefficients into
% the basis coefficients.
a = diag(diag(u))*invu(u, invl(l, p*f));
N = size(a,1);

indices = [find(diff(k)); N]; % index markers delinating degree change between rows
degrees = k(indices);
c = zeros([dim(D,k(end)) size(f,2)]);

row_id = 1;
v_position = 1;
v = v.';
connmat = opt.connection(dim(D,k(end)));

% Degree q space
degree_index = 0;
for q = 0:k(end);
  % If no degree-q information is here, skip:
  if isempty(find(k==q))
    continue
  else
    degree_index = degree_index + 1;
  end
  % Determine the columns of output coefficients associated with this degree
  current_dim = subdim(D, q); % The current size of k-vectors
  cols = dim(D,q-1)+1;
  cols = cols:(cols+current_dim-1);

  % Determine the rows of G associated with this degree
  %rows = row_id:indices(q+1);
  rows = row_id:indices(degree_index);

  degree_indices = v_position:(v_position+length(rows)*current_dim-1);
  degree_indices = reshape(degree_indices.', [current_dim, length(rows)]).';

  Cmat = connmat(cols,cols);
  
  % First compute inv(u)*v:
  v(degree_indices) = invu(u(rows,rows), v(degree_indices));

  M = Cmat*opt.ip(D,q)/Cmat(1,1)^2; % This assumes that you're normalizing the ip
  c(cols,:) = (a(rows,:).'*(v(degree_indices)*M)).';

  row_id = rows(end)+1;
  v_position = degree_indices(end)+1;
end
