function[c] = least_coeffs(l,u,p,v,k,D,f, varargin)
% tensor_coeffs -- Computes basis coefficients for interpolation
%
% c = tensor_coeffs(l,u,p,v,k,dim,f, [ip=alpha!])
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

persistent invu invl dim subdim ip strict_inputs
if isempty(invu)
  from labtools import strict_inputs
  from labtools.linalg import triu_back_substitute as invu
  from labtools.linalg import tril_forward_substitute as invl

  from speclab.common.tensor import space_dimension as dim
  from speclab.common.tensor import subspace_dimension as subdim

  from speclab.orthopoly.interp import monomial_ip as ip
end

opt = strict_inputs({'ip'}, {ip}, [], varargin{:});

% These are the coefficients for the cardinal interpolants under the de Boor
% inner product. The remainder of the code translates these coefficients into
% the basis coefficients.
a = diag(diag(u))*invu(u, invl(l, p*f));
N = size(a,1);

indices = [find(diff(k)); N]; % index markers delinating degree change between rows
c = zeros([dim(D,k(end)) size(f,2)]);

row_id = 1;
v_position = 1;
v = v.';

for q = 0:k(end);
  % Determine the columns of W associated with this degree
  current_dim = subdim(D, q); % The current size of k-vectors
  cols = dim(D,q-1)+1;
  cols = cols:(cols+current_dim-1);

  % Determine the rows of W associated with this degree
  rows = row_id:indices(q+1);

  degree_indices = v_position:(v_position+length(rows)*current_dim-1);
  degree_indices = reshape(degree_indices.', [current_dim, length(rows)]).';
  
  % First compute inv(u)*v:
  v(degree_indices) = invu(u(rows,rows), v(degree_indices));

  % "Orthodox" way
  %c(cols,:) = (a(rows,:).'*W(rows,cols)).';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % About to do some crazy stuff
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  M = opt.ip(D,q);
  c(cols,:) = (a(rows,:).'*(v(degree_indices)*M)).';
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % I think de Boor pulled a fast one on me: the 'orthodox' computation
  % generates coefficients for the scaled multimonomials alpha!*x^alpha. Here
  % I've included the scaling alpha! in the coefficients to make things
  % compatible with mulitmonomial
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  row_id = rows(end)+1;
  v_position = degree_indices(end)+1;
end
