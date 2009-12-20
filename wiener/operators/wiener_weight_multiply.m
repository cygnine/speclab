function[F] = wiener_weight_multiply(F,varargin)
% [F] = wiener_weight_multiply(F, {s=1, t=0, scale=1})
%
%      The input F is a collection of *unweighted* Wiener (s=1,t=0) modal
%      coefficients representing a function f. The output is the modal
%      coefficients vector representing the function f*w^s, where w is the
%      Wiener s=1 phase-shifted square root weight factor. See
%      speclab.wiener.weights.phase_shifted_square_root.

persistent defaults 
if isempty(defaults)
  from speclab.wiener import defaults
end

opt = defaults(varargin{:});

N = length(F);
fconnection = spdiags(ones([N,2]), [0,1], N,N);
for scount = 1:opt.s
  F = fconnection*F;
end

F = F/(sqrt(opt.scale)*(sqrt(2)/i)^opt.s);
