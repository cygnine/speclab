function[dwPsi] = dweighted_fseries(theta,k,varargin)
% [dwpsi] = weighted_fseries(theta,k,{gamma=0, delta=0, shift=0, scale=1,normalization='normal'})
%
%     Evaluates the derivative of the weighted generalized Fourier Series
%     functions with indices k at the locations theta over [0,2*pi]. The matrix
%     psi returned has dimensions length(theta) x length(k).
%
%     The parameters gamma and delta specify the generalization class of the
%     Fourier Series. gamma=delta=0 is the canonical basis. shift and scale
%     define the affine mapping for evaluation of the functions over general
%     intervals. normalization denotes the type of normalization for the
%     trigonometric functions. The default 'normal' indicates that the functions
%     will be L^2-normalized with respect to the unweighted L^2 norm over the
%     interval.

persistent fseries dfseries sqrt_weight dsqrt_weight defaults
if isempty(fseries)
  from speclab.fourier.eval import fseries dfseries
  from speclab.fourier.weights import phase_shifted_sqrt_weight as sqrt_weight
  from speclab.fourier.weights import dphase_shifted_sqrt_weight as dsqrt_weight
  from speclab.fourier import defaults
end

opt = defaults(varargin{:});

theta = theta(:);
N_theta = length(theta);

Psi = fseries(theta,k,opt);
dPsi = dfseries(theta,k,opt);

w = sqrt_weight(theta,opt);
dw = dsqrt_weight(theta,opt);

dwPsi = spdiags(dw,0,N_theta,N_theta)*Psi + ...
       spdiags(w,0,N_theta,N_theta)*dPsi;
