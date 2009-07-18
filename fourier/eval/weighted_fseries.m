function[wPsi] = weighted_fseries(theta,k,varargin)
% [WPSI] = WEIGHTED_FSERIES(THETA,K,{GAMMA=0, DELTA=0, SHIFT=0, SCALE=1,NORMALIZATION='normal'})
%
%     Evaluates the weighted generalized Fourier Series functions with indices K
%     at the locations THETA over [0,2*pi]. The matrix PSI returned has
%     dimensions length(THETA) x length(K).
%
%     The parameters GAMMA and DELTA specify the generalization class of the
%     Fourier Series. GAMMA=DELTA=0 is the canonical basis. SHIFT and SCALE
%     define the affine mapping for evaluation of the functions over general
%     intervals. NORMALIZATION denotes the type of normalization for the
%     trigonometric functions. The default 'normal' indicates that the functions
%     will be L^2-normalized with respect to the unweighted L^2 norm over the
%     interval.

global handles;
fseries = handles.speclab.fourier.eval.fseries;
sqrt_weight = handles.speclab.fourier.weights.phase_shifted_sqrt_weight;
opt = handles.speclab.fourier.defaults(varargin{:});

theta = theta(:);
N_theta = length(theta);

Psi = fseries(theta,k,opt);

w = sqrt_weight(theta,opt);

wPsi = spdiags(w,0,N_theta,N_theta)*Psi;
