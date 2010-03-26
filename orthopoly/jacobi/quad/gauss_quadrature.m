function[x,w] = gauss_quadrature(N,varargin)
% [x,w] = gauss_quadrature(N,{alpha=-1/2,beta=-1/2,shift=0,scale=1})
%
%     Returns the N-point Gaussian quadrature rule for the Jacobi polynomials. 
%     The weight function is given by speclab.orthopoly.jacobi.weights.weight.
%
%     For dimensions > 1, this function returns a *tensorized* grid with N nodes
%     on each dimension.

persistent jac pss gq tensorize tensorize_vectors npre npost defaults recurrence
if isempty(jac)
  from speclab.orthopoly.jacobi import defaults
  from speclab.common import physical_scaleshift as pss
  from speclab.common.tensor import tensorize_vectors tensorize
  from speclab.orthopoly import gauss_quadrature as gq
  from speclab.orthopoly.jacobi.coefficients import recurrence
end

opt = defaults(varargin{:});
[alpha,beta,scale,shift] = deal(opt.alpha,opt.beta,opt.scale,opt.shift);

tol = 1e-8;

% If all bases are Chebyshev, we can do this a simple way:
if all(abs(alpha+1/2)<tol) && all(abs(beta+1/2)<tol);
  w = pi/N*ones([N 1]);
  temp = linspace(pi,0,N+1)';
  x = cos(temp(1:N)-pi/(2*N));

  if opt.dim>1  % copy+paste from speclab.orthopoly.gauss_quadrature
    x = tensorize(x, opt.dim);
    w = tensorize(w, opt.dim);
    w = prod(w, 2);
  end

% If all bases are the same, only call gq once
elseif all(abs(alpha-alpha(1))<tol) && all(abs(beta-beta(1))<tol)
  [a,b] = recurrence(N+1,'alpha',opt.alpha(1), 'beta', opt.beta(1));

  [x,w] = gq(a,b,N,'dim', opt.dim);

% *sigh*, everything is sufficiently different
else
  xs = cell([opt.dim 1]);
  ws = cell([opt.dim 1]);

  for q = 1:opt.dim
    [a,b] = recurrence(N+1, 'alpha', opt.alpha(q), 'beta', opt.beta(q));

    [xs{q}, ws{q}] = gq(a,b,N);
  end

  x = tensorize_vectors(xs{:});
  w = prod(tensorize_vectors(ws{:}), 2);
end

switch lower(opt.weight_normalization)
case 'probability'
  for q = 1:opt.dim
    [a,b] = recurrence(1,'alpha',opt.alpha(q),'beta',opt.beta(q));
    w = w/b;
  end
case {'lebesgue', 'lebesque'}
  for q = 1:opt.dim
    if (opt.alpha(q)==0) & (opt.beta(q)==0)
      w = opt.scale(q)/2;
    else
      error('''Lebesgue'' weight normalization not supported for non-Legendre expansions');
    end
  end
otherwise
  %error('Unrecognized weight specification');
end

x = pss(x,opt);
