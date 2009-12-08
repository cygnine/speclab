function[p] = eval_jacobi_poly(x,n,varargin)
% [p] = eval_jacobi_poly(x,n,{alpha=-1/2,beta=-1/2,d=0,dim=1,shift=0,scale=1})
%
%     Evaluates the nth generalized Jacobi polynomials at the locations x. 
%     The weight function is (1-1/scale*(x-shift))^alpha*(1+1/scale*(x-shift))^beta
%     Is vectorized in x and n.
%
%     The optional dimension dim specifies the dimension of the x-variable. If
%     dim > 1, a tensor-product function is assumed. If dim > 1, x is assumed to
%     be an M x dim array, with each row corresponding to a d-dimensional data
%     point. In the multidimensional case the linear index n is given by the
%     array indexing procedure from
%     speclab.common.tensor.linear_to_array_indexing.

persistent recurrence defaults eval_polynomial indexing npost npre
if isempty(recurrence)
  from speclab.orthopoly1d.jacobi.coefficients import recurrence
  from speclab.orthopoly1d.jacobi import defaults
  from speclab.orthopoly1d import eval_polynomial
  from speclab.common.tensor import linear_to_array_indexing as indexing

  from speclab.common.tensor import node_preprocessing as npre
end

opt = defaults(varargin{:});

[x,info] = npre(x,opt.dim);

if opt.dim==1
  n_array = n(:); % no need to call linear_to_array_indexing
  N = max(n)+2;
else
  n_array = indexing(n,'dim',opt.dim);
  N = max(max(n_array)) + 2;
end


switch opt.normalization
case 'normal'
  % Loop over each dimension, although if memory weren't a concern, we could do
  % this all in one step
  one_d_opts.d = opt.d;
  one_d_opts.normalization = opt.normalization;
  p = ones([size(x,1), size(n_array,1)]);

  factor = 1;

  for q = 1:opt.dim
    [a,b] = recurrence(N+1,'alpha', opt.alpha(q), 'beta', opt.beta(q));

    factor = factor*b(1);

    one_d_opts.shift = opt.shift(q);
    one_d_opts.scale = opt.scale(q);
    p = p.*eval_polynomial(x(:,q),a,b,n_array(:,q),one_d_opts);
  end

  switch opt.weight_normalization
  case ''
    % do nothing
  case 'probability'
    p = p*sqrt(factor);
  otherwise
    error('Weight normalization type not supported');
  end

case 'monic'
  error('Not yet implemented');
otherwise
  error('Normalization type not supported');
end
