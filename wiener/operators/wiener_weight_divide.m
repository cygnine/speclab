function[F] = wiener_weight_divide(F,varargin)
% [F] = wiener_weight_divide(F, {s=1, t=0, scale=1})
%
%      The input F is a collection of *unweighted* Wiener (s=1,t=0) modal
%      coefficients representing a function f. The output is the modal
%      coefficients vector representing the function f/w^s, where w is the
%      Wiener s=1 phase-shifted square root weight factor. See
%      speclab.wiener.weights.phase_shifted_square_root.
%
%      Since I can only see how to do this with matlab for-loops, it's amazingly
%      slow, even though it would be speedy in a compiled language. The
%      bottleneck is really triu_sparse_invert.
%
%      TODO: optimize via sliding s-length updater over vector (removes outer
%      for loop).

persistent defaults matinv
if isempty(defaults)
  from speclab.wiener import defaults
  from labtools.linalg import triu_sparse_invert as matinv
end

opt = defaults(varargin{:});

N = length(F);
fconnection = spdiags(ones([N,2]), [0,1], N,N);
for scount = 1:opt.s
  F = matinv(fconnection, F, 'bandwidth', 2);
end

F = F*sqrt(opt.scale)*(sqrt(2)/i)^opt.s;
