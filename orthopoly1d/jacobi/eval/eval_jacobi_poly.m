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

persistent recurrence defaults driver
if isempty(recurrence)
  from speclab.orthopoly1d.jacobi.coefficients import recurrence
  from speclab.orthopoly1d.jacobi import defaults
  from speclab.orthopoly1d import eval_poly_driver as driver
end

opt = defaults(varargin{:});

poly_parameters = struct('alpha', opt.alpha, 'beta', opt.beta);

p = driver(x,n,opt.d,recurrence,opt.dim,opt.shift,opt.scale,opt.normalization, ...
           opt.weight_normalization, poly_parameters);
