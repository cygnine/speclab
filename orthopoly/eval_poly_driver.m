function[p] = eval_poly_driver(x,n,d,recurrence,dim,shift,scale,normalization,weight_normalization,poly_parameters,associated_index)
% [p] = eval_poly_driver(x,n,d,recurrence,dim,shift,scale,normalization,weight_normaliztion,poly_parameters, associated_index)
%
%     Evaluates the nth degree orthogonal polynomials at the locations x.  Is
%     vectorized in x and n. 
%
%     The derivative counter d specifies which derivatives to take and the size
%     of this input vector corresponds to the length of the third dimension of
%     the output p. The output p is of size length(x) x length(n) x length(d).
%
%     The input 'recurrence' is a function handle pointing to the recurrence
%     relation for the orthogonal polynomial family.
%
%     The dimension 'dim' specifies the dimension of the x-variable. If
%     dim > 1, a tensor-product function is assumed. If dim > 1, x is assumed to
%     be an M x dim array, with each row corresponding to a d-dimensional data
%     point. In the multidimensional case the linear index n is given by the
%     array indexing procedure from
%     speclab.common.tensor.linear_to_array_indexing.
%
%     The 'shift' and 'scale' inputs are the multidimensional affine shift and
%     scale parameters relative to the standard interval.
%
%     'poly_parameters' is a struct with fields given by any necessary
%     parameters for the orthogonal polynomial family.
%
%     'normalization' and 'weight_normalization' are strings that determine how
%     the polynomials are normalized ('normalization') with respect to the
%     particular form of the weight_function ('weight_normalization').
%
%     'associated_index' takes values 0, 1, 2, ... and indicates the
%     `associated' polynomials to be evaluated. With value 0, they are the usual
%     orthogonal polynomials. They are generated by shifting the recurrence
%     coefficients by the indicated number of values.
%
%     This is the main preprocessing driver for all tensor-product orthogonal
%     polynomial routines.

persistent eval_polynomial indexing npre
if isempty(eval_polynomial)
  from speclab.orthopoly import eval_polynomial
  from speclab.common.tensor import linear_to_array_indexing as indexing
  from speclab.common.tensor import node_preprocessing as npre
end

[x,info] = npre(x,dim);

if dim==1
  n_array = n(:); % no need to call linear_to_array_indexing
  N = max(n)+2;
else
  n_array = indexing(n,'dim',dim);
  N = max(max(n_array)) + 2;
end

% Loop over each dimension, although if memory weren't a concern, we could do
% this all in one step
one_d_opts.d = d;
one_d_opts.normalization = normalization;
p = ones([size(x,1), size(n_array,1)]);
factor = 1; % Used for the `probability' weight normalization

for q = 1:dim
  [a,b] = recurrence(N+associated_index+1,poly_parameters(q));
  a(1:associated_index) = [];
  b(1:associated_index) = [];

  factor = factor*b(1);

  one_d_opts.shift = shift(q);
  one_d_opts.scale = scale(q);
  p = p.*eval_polynomial(x(:,q),a,b,n_array(:,q),one_d_opts);
end

% Re-normalize if strange weight is used
switch weight_normalization
case ''
  % do nothing
case 'probability'
  p = p*sqrt(factor);
otherwise
  error('Weight normalization type not supported');
end
