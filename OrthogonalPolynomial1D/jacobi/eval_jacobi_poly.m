function[p] = eval_jacobi_poly(x,n,varargin)
% [P] = EVAL_JACOBI_POLY(X,N,{ALPHA=-1/2,BETA=-1/2,SHIFT=0,SCALE=1})
%
%     Evaluates the nth generalized Jacobi polynomials at the locations x. 
%     The weight function is (1-1/scale*(x-shift))^alpha*(1+1/scale*(x-shift))^beta
%     Is vectorized in x and n.
%
% 20080524: acn

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
sss = handles.speclab.common.standard_scaleshift_1d;
jac.defaults;

N = max(n)+2;

[a,b] = jacobi_recurrence(N+1,opt);
x = sss(x,opt.scale,opt,shift);

switch opt.normalization
case 'normal'
  p = eval_normalized_orthogonal_poly(x,a,b,n,opt);
  p = p/sqrt(opt.scale);
end
