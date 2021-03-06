function[wPsi] = weighted_fseries(theta,k,varargin)
% [wpsi] = weighted_fseries(theta,k,{gamma=0, delta=0, shift=0, scale=1,normalization='normal'})
%
%     Evaluates the weighted generalized Fourier Series functions with indices K
%     at the locations theta over [-pi,pi]. The matrix PSI returned has
%     dimensions length(theta) x length(k).
%
%     The parameters gamma and delta specify the generalization class of the
%     Fourier Series. gamma=delta=0 is the canonical basis. shift and scale
%     define the affine mapping for evaluation of the functions over general
%     intervals. normalization denotes the type of normalization for the
%     trigonometric functions. The default 'normal' indicates that the functions
%     will be L^2-normalized with respect to the unweighted L^2 norm over the
%     interval.

persistent fseries sqrt_weight defaults
if isempty(fseries)
  from speclab.fourier.eval import fseries
  from speclab.fourie.weights import phase_shifted_sqrt_weight as sqrt_weight
  from speclab.fourier import defaults
end

opt = defaults(varargin{:});

theta = theta(:);
N_theta = length(theta);

Psi = fseries(theta,k,opt);

w = sqrt_weight(theta,opt);

wPsi = spdiags(w,0,N_theta,N_theta)*Psi;
