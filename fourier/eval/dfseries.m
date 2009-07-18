function[dPsi] = dfseries(theta,k,varargin)
% [DPSI] = DFSERIES(THETA,K,{GAMMA=0, DELTA=0, SHIFT=0, SCALE=1,NORMALIZATION='normal'})
%
%     Evaluates the derivative of the generalized Fourier Series functions with
%     indices K at the locations THETA over [0,2*pi]. The matrix DPSI returned has
%     dimensions length(THETA) x length(K).
%
%     The parameters GAMMA and DELTA specify the generalization class of the
%     Fourier Series. GAMMA=DELTA=0 is the canonical basis. SHIFT and SCALE
%     define the affine mapping for evaluation of the functions over general
%     intervals. NORMALIZATION denotes the type of normalization for the
%     trigonometric functions (note that the normalization specifications refers
%     to the functions themselves, not the derivatives). The default 'normal'
%     indicates that the functions will be L^2-normalized, where L^2 refers to
%     the weighted norm specified by GAMMA and DELTA. See the files in
%     speclab/fourier/weights.

global handles;
jac = handles.speclab.OrthogonalPolynomial1D.jacobi;
rtheta = handles.speclab.fourier.maps;
opt = handles.speclab.fourier.defaults(varargin{:});

theta = theta(:);
N_theta = length(theta);
k = k(:);
N_k = length(k);

k_not_0 = (k~=0);
k_is_0 = ~(k_not_0);

alpha = opt.delta - 1/2;
beta = opt.gamma - 1/2;

r = rtheta.theta_to_r(theta,opt);
rc = rtheta.theta_to_rc(theta,opt);
dr = rtheta.dr_dtheta(theta,opt);
drc = rtheta.drc_dtheta(theta,opt);

dPsi = zeros([N_theta, N_k]);

dPsi = 1/2*spdiags(dr,0,N_theta,N_theta)*...
       jac.eval_jacobi_poly(r,abs(k),'alpha', alpha, ...
                                     'beta',  beta, ...
                                     'd', 1);
dPsi(:,k_is_0) = dPsi(:,k_is_0)*sqrt(2);

Psi(:,k_is_0) = 1/sqrt(2)*p1(:,k_is_0);

if any(k_not_0)
  p2 = jac.eval_jacobi_poly(r,abs(k(k_not_0)), 'alpha', alpha+1, ...
                                               'beta',  beta+1);
  dp2 = jac.eval_jacobi_poly(r,abs(k(k_not_0)), 'alpha', alpha+1, ...
                                                'beta',  beta+1, ...
                                                'd', 1);
  odd_term = i*spdiags(drc,0,N_theta,N_theta)*p2 + ...
             i*spdiags(rc.*dr,0,N_theta,N_theta)*dp2;
  dPsi(:,k_not_0) = dPsi(:,k_not_0) + ...
                    1/2*odd_term*spdiags(sign(k(k_not_0)),0,N_k,N_k);
end