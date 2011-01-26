function[x,w] = gauss_quadrature(N,varargin)
% gauss_quadrature -- Gauss quadrature rule for Hermite polynomials
%
% [x,w] = gauss_quadrature(N,{mu=0,shift=0,scale=1})
%
%     Returns the N-point Gaussian quadrature rule for the Hermite polynomials. 
%     The weight function is given by speclab.orthopoly.hermite.weights.weight.

persistent recurrence gq pss defaults tensorize tensorize_vectors
if isempty(recurrence)
  from speclab.orthopoly.hermite import defaults recurrence
  from speclab.orthopoly import gauss_quadrature as gq
  from speclab.common import physical_scaleshift as pss
  from speclab.common.tensor import tensorize_vectors tensorize
end

opt = defaults(varargin{:});

xs = cell([opt.dim 1]);
ws = cell([opt.dim 1]);

for q = 1:opt.dim
  [a,b] = recurrence(N+1, 'mu', opt.mu(q));

  [xs{q}, ws{q}] = gq(a,b,N);
end

x = tensorize_vectors(xs{:});
w = prod(tensorize_vectors(ws{:}), 2);

%[a,b] = recurrence(N+1,opt);
%[x,w] = gq(a,b,N);

switch opt.weight_normalization
case 'probability'
  for q = 1:opt.dim
    [a,b] = recurrence(1,'mu',opt.mu(q));
    w = w/b;
  end
  %w = w/b(1);
end

x = pss(x,opt);
