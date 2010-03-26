function[c] = opoly_li_coeffs(l,W,p,k,u,D,f)
% opoly_li_coeffs -- Computes orthogonal polynomial coefficients for interpolation
%
% c = opoly_li_coeffs(l,W,p,k,u,dim,f)
%
%     Using the 5-tuple (l,W,p,k,u) that is output from de Boor's algorithm
%     (modified by J. Li) opoly_li_lu, this function uses the point-evaluations
%     f to compute the orthogonal polynomial coefficients necessary for
%     interpolation.  Obviously, f must be an N-vector, where N is the number of
%     rows of W (or l, or p, or u).  Unfortunately, the dimension of the problem
%     cannot be divined from the 5-tuple, so the extra input dim denoting the
%     dimension is necessary.
%
%     The interpolant may be evaluated as 
%
%       p(x) = sum_n c(n) * m(x,n),
%
%     where m(x,n) is a function that evaluates the n'th orthogonal polynomial
%     at the location x. (e.g. speclab.orthopoly.jacobi.eval.eval_jacobi_poly
%     is a possibility). The output c has length size(W,2). Note that due to J.
%     Li's choice of inner product, the particular class of orthogonal
%     polynomials is irrelevant.
%
%     Note that this function is vectorized in the columns of f, meaning that if
%     f is a collection of column vectors, this function returns a collection of
%     column vectors containing the monomial representation.

persistent invu invl dim subdim ip
if isempty(invu)
  from labtools.linalg import triu_back_substitute as invu
  from labtools.linalg import tril_forward_substitute as invl

  from speclab.common.tensor import space_dimension as dim
  from speclab.common.tensor import subspace_dimension as subdim

  ip = @(d,k) speye(subdim(d,k));
  %from speclab.orthopoly.interp import monomial_ip as ip
end

% These are the coefficients for the cardinal interpolants under the de Boor
% inner product. The remainder of the code translates these coefficients into
% multimonomial coefficients.
a = diag(diag(u))*invu(u, invl(l, p*f));
N = size(a,1);

W = invu(u, W);

c = zeros([size(W,2) size(f,2)]);
kmax = max(k);
kdiff = diff(k);
indices = [find(kdiff); N];

% Deal with k==0 separately:
c(1,:) = a(1,:)*W(1,1);
row_id = 2;

% For all k>0:
for q = 1:kmax;
  % Determine the columns of W associated with this degree
  current_dim = subdim(D, q); % The current size of k-vectors
  cols = dim(D,q-1)+1;
  cols = cols:(cols+current_dim-1);

  % Determine the rows of W associated with this degree
  rows = row_id:indices(q+1);

  % "Orthodox" way
  %c(cols,:) = (a(rows,:).'*W(rows,cols)).';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % About to do some crazy stuff
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  M = ip(D, q);
  c(cols,:) = (a(rows,:).'*(W(rows,cols)*M)).';
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % I think de Boor pulled a fast one on me: the 'orthodox' computation
  % generates coefficients for the scaled multimonomials alpha!*x^alpha. Here
  % I've included the scaling alpha! in the coefficients to make things
  % compatible with mulitmonomial
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  row_id = rows(end)+1;
end
