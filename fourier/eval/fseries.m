function[Psi] = fseries(theta,k,varargin)
% [PSI] = FSERIES(THETA,K,{GAMMA=0, DELTA=0, SHIFT=0, SCALE=1,NORMALIZATION='normal'})
%
%     Evaluates the generalized Fourier Series functions with indices K at the
%     locations THETA over [0,2*pi]. The matrix PSI returned has dimensions
%     length(THETA) x length(K).
%
%     The parameters GAMMA and DELTA specify the generalization class of the
%     Fourier Series. GAMMA=DELTA=0 is the canonical basis. SHIFT and SCALE
%     define the affine mapping for evaluation of the functions over general
%     intervals. NORMALIZATION denotes the type of normalization for the
%     trigonometric functions. The default 'normal' indicates that the functions
%     will be L^2-normalized, where L^2 refers to the weighted norm specified by
%     GAMMA and DELTA. See the files in speclab/fourier/weights.

global handles;
jac = handles.speclab.orthopoly1d.jacobi;
fourier = handles.speclab.fourier;
rtheta = handles.speclab.fourier.maps;
opt = handles.speclab.fourier.defaults(varargin{:});

theta = theta(:);
N_theta = length(theta);
k = k(:);
N_k = length(k);

if fourier.classic_fourier(opt);
  sss = handles.speclab.common.standard_scaleshift_1d;
  Psi = 1/sqrt(2*pi)*exp(i*sss(theta,opt)*k.');
  return;
end

k_not_0 = (k~=0);
k_is_0 = ~(k_not_0);

N_k_not_0 = N_k - sum(k_is_0);

alpha = opt.delta - 1/2;
beta = opt.gamma - 1/2;

r = rtheta.theta_to_r(theta,opt);
rc = rtheta.theta_to_rc(theta,opt);

Psi = zeros([N_theta, N_k]);

p1 = jac.eval.eval_jacobi_poly(r,abs(k),'alpha', alpha, ...
                                        'beta',  beta);

Psi(:,k_is_0) = 1/sqrt(2)*p1(:,k_is_0);

if any(k_not_0)
  p2 = jac.eval.eval_jacobi_poly(r,abs(k(k_not_0))-1, 'alpha', alpha+1, ...
                                                      'beta',  beta+1);
  p2 = i*spdiags(rc,0,N_theta,N_theta)*p2...
       *spdiags(sign(k(k_not_0)), 0,N_k_not_0, N_k_not_0);
  
  Psi(:,k_not_0) = 1/2*(p1(:,k_not_0) + p2);
end
