function[modes] = piecewise_poly_modes(monomial_coefficients,mesh,N);
% [MODES] = PIECEWISE_POLY_MODES(MONOMIAL_COEFFICIENTS,MESH,N);
%
%     Compute exact Fourier coefficients given a piecewise polynomials. MESH is
%     a vector of length k; MONOMIAL_COEFFICIENTS is an N x k matrix of modal
%     coefficients where the local basis on each element is unnormalized
%     monomials of degree N-1 or less. Assumes the interval of integration is
%     [0,2*pi]. This function returns the first N Fourier modal coefficients
%     exactly computed from the piecewise polynomial representation. The indices
%     of the N modal coefficeints are determined by canonical modal numbering
%     scheme. See integer_range.

global handles;
fourier = handles.speclab.fourier;

ks = fourier.integer_range(N);
[N,k] = size(monomial_coefficients);
