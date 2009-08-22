function[modes] = matlab_modes_pre_ifft(modes);
% [modes] = matlab_modes_pre_ifft(modes);
%
%     The inverse of matlab_fft_to_modes. Prepares the modal coefficients
%     corresponding to the basis expansion 1/sqrt(2*pi)*exp(i*k*theta) so that
%     they can be input into matlab's ifft to give back nodal evaluations at
%     equidistant points.

modes = ifftshift(modes)*length(modes)/sqrt(2*pi);
