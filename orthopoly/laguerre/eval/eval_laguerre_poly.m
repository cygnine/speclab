function[p] = eval_laguerre_poly(x,n,varargin)
% eval_laguerre_poly -- evaluate Laguerre polynomials
%
% [p] = eval_laguerre_poly(x,n,{alpha=0,d=0,shift=0,scale=1})
%
%     Evaluates the nth generalized Laguerre polynomials at the locations x.
%     This function is vectorized in x and n. The interval of approximation is
%     scaled to interval*scale + shift. The Jacobian resulting from this affine
%     transform is built into the weight function.
%
%     The optional dimension dim specifies the dimension of the x-variable. If
%     dim > 1, a tensor-product function is assumed. If dim > 1, x is assumed to
%     be an M x dim array, with each row corresponding to a d-dimensional data
%     point. In the multidimensional case the linear index n is given by the
%     array indexing procedure from
%     speclab.common.tensor.linear_to_array_indexing.

persistent recurrence defaults driver
if isempty(recurrence)
  from speclab.orthopoly.laguerre.coefficients import recurrence
  from speclab.orthopoly.laguerre import defaults
  from speclab.orthopoly import eval_poly_driver as driver
end

opt = defaults(varargin{:});

poly_parameters = struct('alpha', opt.alpha);

p = driver(x,n,opt.d,recurrence,opt.dim,opt.shift,opt.scale,opt.normalization, ...
           opt.weight_normalization, poly_parameters);
