function[modes] = matlab_modes_pre_ifft(modes);
% [MODES] = MATLAB_MODES_PRE_IFFT(MODES);
%
%     The inverse of MATLAB_FFT_TO_MODES. Prepares the modal coefficients
%     corresponding to the basis expansion 1/sqrt(2*pi)*exp(i*k*theta) so that
%     they can be input into matlab's ifft to give back nodal evaluations at
%     equidistant points.

modes = ifftshift(modes)*length(modes)/sqrt(2*pi);
