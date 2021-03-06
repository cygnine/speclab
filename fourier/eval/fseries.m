function[Psi] = fseries(theta,k,varargin)
% [psi] = fseries(theta,k,{gamma=0, delta=0, shift=0, scale=1,normalization='normal'})
%
%     Evaluates the generalized Fourier Series functions with indices K at the
%     locations theta over [0,2*pi]. The matrix psi returned has dimensions
%     length(theta) x length(k).
%
%     The parameters gamma and delta specify the generalization class of the
%     Fourier Series. gamma=delta=0 is the canonical basis. shift and scale
%     define the affine mapping for evaluation of the functions over general
%     intervals. normalization denotes the type of normalization for the
%     trigonometric functions. The default 'normal' indicates that the functions
%     will be L^2-normalized, where L^2 refers to the weighted norm specified by
%     gamma and delta. See the files in speclab/fourier/weights.

persistent eval_jacobi_poly classic_fourier theta_to_r theta_to_rc sss defaults
if isempty(eval_jacobi_poly)
  from speclab.orthopoly.jacobi.eval import eval_jacobi_poly
  from speclab.fourier import classic_fourier defaults
  from speclab.fourier.maps import theta_to_r theta_to_rc
  from speclab.common import standard_scaleshift_1d as sss
end

opt = defaults(varargin{:});

theta = theta(:);
N_theta = length(theta);
k = k(:);
N_k = length(k);

if classic_fourier(opt);
  Psi = 1/sqrt(2*pi)*exp(i*sss(theta,opt)*k.');
  return;
end

k_not_0 = (k~=0);
k_is_0 = ~(k_not_0);

N_k_not_0 = N_k - sum(k_is_0);

alpha = opt.delta - 1/2;
beta = opt.gamma - 1/2;

r = theta_to_r(theta,opt);
rc = theta_to_rc(theta,opt);

Psi = zeros([N_theta, N_k]);

p1 = eval_jacobi_poly(r,abs(k),'alpha', alpha, ...
                                        'beta',  beta);

Psi(:,k_is_0) = 1/sqrt(2)*p1(:,k_is_0);

if any(k_not_0)
  p2 = eval_jacobi_poly(r,abs(k(k_not_0))-1, 'alpha', alpha+1, ...
                                                      'beta',  beta+1);
  p2 = i*spdiags(rc,0,N_theta,N_theta)*p2...
       *spdiags(sign(k(k_not_0)), 0,N_k_not_0, N_k_not_0);
  
  Psi(:,k_not_0) = 1/2*(p1(:,k_not_0) + p2);
end
