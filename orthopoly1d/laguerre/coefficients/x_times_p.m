function[mu] = x_times_p(n,alpha,varargin)
% [mu] = x_times_p(n,alpha,{normalization='normal'})
%
%     Computes the connection coefficients between x x L and L. I.e., given
%     the vector n of whole-number indices, and valid Laguerre class parameter
%     alpha, mu is a length(n) x 2 matrix, where each row contains the
%     coefficients defined by the relation:
%
%     x*L^(alpha)_n = mu(1)*L^(alpha-1)_n +
%                     mu(2)*L^(alpha-1)_{n+1}
%
%     Note that alpha must be larger than 0.
%     The optional input normalization defines the normalization of the
%     polynomials L, which affects the values of the parameters mu.

global packages;
opt = packages.speclab.orthopoly1d.laguerre.defaults(varargin{:});
%opt = packages.common.input_schema({'normalization','scale'}, {'normal',1}, [],varargin{:});
n = n(:);
N = length(n);

mu = zeros([N,2]);

if strcmpi(opt.normalization,'normal')
  mu(:,1) = sqrt(n+alpha);
  mu(:,2) = sqrt(n+1);
  %mu = mu/opt.scale;
  mu = mu*opt.scale;
elseif strcmpi(opt.normalization, 'classical')
  mu(:,1) = n+alpha;
  mu(:,2) = -n;
  mu = mu*opt.scale;
else
  fprintf('Error: normalization %s not supported', opt.normalization);
  mu = false;
  return
end
