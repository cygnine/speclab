function[w] = weight(r,varargin)
% [w] = weight(r,{alpha=-1/2,beta=-1/2,shift=0,scale=1,weight_normalization=[]})
%
%     Evaluates the Jacobi polynomial weight function for class (alpha,beta)
%     under the affine map r -> r*scale + shift. When scaled, this weight
%     function absorbs the Jacobian, leaving the point evaluations of the
%     polynomials unchanged with respect to the input 'scale'.

persistent defaults sss
if isempty(defaults)
  from speclab.orthopoly1d.jacobi import defaults
  from speclab.common import standard_scaleshift_1d as sss
end

opt = defaults(varargin{:});
r = sss(r,opt);

%for q = 1:opt.dim
w = (1-r).^opt.alpha;
w = w.*(1+r).^opt.beta;

switch opt.weight_normalization
case 'probability'
  recur = from_as('speclab.orthopoly1d.jacobi.coefficients', 'recurrence');
  [a,b] = recur(1,'alpha',opt.alpha,'beta',opt.beta);
  w = w/(b*opt.scale);
otherwise
  w = w/opt.scale;
end
