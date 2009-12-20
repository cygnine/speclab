function[dPsi] = dfseries(theta,k,varargin)
% [dpsi] = dfseries(theta,k,{gamma=0, delta=0, shift=0, scale=1,normalization='normal'})
%
%     Evaluates the derivative of the generalized Fourier Series functions with
%     indices k at the locations theta over [-pi,pi]. The matrix DPSI returned has
%     dimensions length(theta) x length(k).
%
%     The parameters gamma and delta specify the generalization class of the
%     Fourier Series. gamma=delta=0 is the canonical basis. shift and scale
%     define the affine mapping for evaluation of the functions over general
%     intervals. normalization denotes the type of normalization for the
%     trigonometric functions (note that the normalization specifications refers
%     to the functions themselves, not the derivatives). The default 'normal'
%     indicates that the functions will be L^2-normalized, where L^2 refers to
%     the weighted norm specified by gamma and delta. See the files in
%     speclab/fourier/weights.

persistent rtheta eval_jacobi_poly defaults
if isempty(eval_jacobi_poly)
  from speclab.fourier import maps as rtheta
  from speclab.orthopoly1d.jacobi.eval import eval_jacobi_poly
  from speclab.fourier import defaults
end

opt = defaults(varargin{:});

theta = theta(:);
N_theta = length(theta);
k = k(:);
N_k = length(k);

k_not_0 = (k~=0);
k_is_0 = ~(k_not_0);

N_k_not_0 = N_k - sum(k_is_0);

alpha = opt.delta - 1/2;
beta = opt.gamma - 1/2;

r = rtheta.theta_to_r(theta,opt);
rc = rtheta.theta_to_rc(theta,opt);
dr = rtheta.dr_dtheta(theta,opt);
drc = rtheta.drc_dtheta(theta,opt);

dPsi = zeros([N_theta, N_k]);

dPsi = 1/2*spdiags(dr,0,N_theta,N_theta)*...
       eval_jacobi_poly(r,abs(k),'alpha', alpha, ...
                                          'beta',  beta, ...
                                          'd', 1);
dPsi(:,k_is_0) = dPsi(:,k_is_0)*sqrt(2);

if any(k_not_0)
  p2 = eval_jacobi_poly(r,abs(k(k_not_0))-1, 'alpha', alpha+1, ...
                                                      'beta',  beta+1);
  dp2 = eval_jacobi_poly(r,abs(k(k_not_0))-1, 'alpha', alpha+1, ...
                                                       'beta',  beta+1, ...
                                                       'd', 1);
  odd_term = i*spdiags(drc,0,N_theta,N_theta)*p2 + ...
             i*spdiags(rc.*dr,0,N_theta,N_theta)*dp2;
  dPsi(:,k_not_0) = dPsi(:,k_not_0) + ...
                    1/2*odd_term*spdiags(sign(k(k_not_0)),0,N_k_not_0,N_k_not_0);
end
