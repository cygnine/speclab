function[df] = derivative(f, varargin)
% derivative -- Calculates a Fourier derivative using the FFT
%
% df = derivative(f, {alpha=1,scale=1,shift=0})
%
%     Computes the deriative of the pointwise values f (at the same points)
%     using two FFT's. The nodal locations are assumed to be the canonical
%     Fourier locations. 
%
%     alpha is the order of the derivative. ANY real number is supported, using
%     the Fourier Series definition of the derivative (or interpolation spaces).
%
%     The operation is performed along the first dimension -- other dimensions
%     are ignored.
%
%     You can write an optimized version of this function for a fixed N by using
%     the (i)ffft_online/overhead routines (also so that you don't have to use
%     strict_inputs).

persistent ffft iffft integer_range spdiag strict_inputs
if isempty(ffft)
  from speclab.fourier.fft import ffft iffft
  from speclab.common import integer_range
  from labtools import spdiag strict_inputs
end

opt = strict_inputs({'alpha', 'shift', 'scale'}, {1, 0, 1}, [], varargin{:});
temp = opt;
temp = rmfield(temp, 'alpha');

N = size(f,1);
fourier_factors = spdiag((i*integer_range(N)).^opt.alpha);

df = iffft(fourier_factors*ffft(f, temp), temp);
