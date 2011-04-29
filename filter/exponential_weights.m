function[factors] = exponential_weights(eta,varargin)
% exponential_weights -- attenuation factors for exponential filter
%
% factors = exponential_weights(eta,{alpha=-log(eps), p=8, preservation_fraction=0.5})
%
%     For 0 <= eta <= 1 (array), computes the attentuation factor for the
%     exponential filter of order p, with terminal attentuation exp(-alpha).
%     Preserves modes with 0 <= eta <= preservation_fraction.

persistent strict_inputs
if isempty(strict_inputs)
  from labtools import strict_inputs
end

inputs = {'alpha', 'p', 'preservation_fraction'};
defaults = {-log(eps), 8, 0.5};
opt = strict_inputs(inputs, defaults, [], varargin{:});

pflags = eta<=opt.preservation_fraction;

factors = ones(size(eta));
factors(~pflags) = exp(-opt.alpha*((eta(~pflags) - opt.preservation_fraction)./...
                                   (1 - opt.preservation_fraction)).^opt.p);
