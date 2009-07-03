function[modes] = matlab_fft_to_modes(modes);
% [MODES] = MATLAB_FFT_TO_MODES(MODES);
%
%     The input MODES is from the direct output of Matlab's FFT. This function
%     scales the modes so that they correspond to modal coefficients in the
%     basis expansion 1/sqrt(2*pi)*exp(i*k*theta), where k is an integer index,
%     negatively biased when N = length(modes) is even.

modes = fftshift(modes)/length(modes)*sqrt(2*pi);
