function[c] = jacobi_to_monomial(c, cell_boundaries, varargin)
% jacobi_to_monomial -- Jacobi polynomial to monomial basis connection
%
% c = jacobi_to_monomial(c, cell_boundaries, {alpha=0, beta=0})
%
%     Translates cellwise local Jacobi polynomial expansion coefficients into
%     modal coefficients for the global monomial basis x^n. The coefficient
%     matrix c is an N x K matrix, representing an N-1 degree polynomial on each
%     of K cells. The coefficients in each column correspond to modal
%     coefficients for L^2-normalized Jacobi polynomials of order (alpha, beta). 
%
%     The length-(K+1) vector cell_boundaries specifies the global vertices of
%     the cells. This is necessary to provide monomial global coefficients. 

persistent input_schema eval_jacobi_poly gq
if isempty(input_schema)
  from labtools import input_schema
  from speclab.orthopoly.jacobi.eval import eval_jacobi_poly
  from speclab.orthopoly.jacobi.quad import gauss_quadrature as gq
end

opt = input_schema({'alpha', 'beta'}, {0,0}, [], varargin{:});

[N,K] = size(c);
% Rudimentary error checking
if length(cell_boundaries)~=(K+1)
  error('The cell_boundaries input must have length==(1 + # of columns of c)');
end

[r,w] = gq(N, opt);
ps = eval_jacobi_poly(r,0:(N-1),opt);

% For each cell: 
%    - generate local nodes
%    - evaluate global monomials at those nodes
%    - form connection matrix
%    - apply connection matrix to appropriate column of c
for q = 1:K;
  cell_shift = mean(cell_boundaries(q:(q+1)));
  cell_scale = (cell_boundaries(q+1) - cell_boundaries(q))/2;

  x = r*cell_scale + cell_shift;

  % Connection matrix
  connection = zeros(N);
  % yeah, yeah, this can be vectorized
  for qq = 1:N
    fx = x.^(qq-1);
    connection(:,qq) = (ps'*(w.*fx));
  end

  c(:,q) = inv(connection)*c(:,q);
end
