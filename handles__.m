function[hs] = handles__()
% [HS] = HANDLES__()
%
%     Returns directory pointers for common module

% This is by default
hs.base = fileparts(mfilename('fullpath'));

hs.common.base = fullfile(hs.base,'common');
hs.fourier.base = fullfile(hs.base,'fourier');
  hs.fourier.eval.base = fullfile(hs.fourier.base, 'eval');
  hs.fourier.quad.base = fullfile(hs.fourier.base, 'eval');
  hs.fourier.weights.base = fullfile(hs.fourier.base, 'eval');
  hs.fourier.maps.base = fullfile(hs.fourier.base, 'eval');
hs.monomials.base = fullfile(hs.base,'monomials');
hs.NewtonPolynomials.base = fullfile(hs.base,'NewtonPolynomials');

hs.OrthogonalPolynomial1D.base = fullfile(hs.base,'OrthogonalPolynomial1D');
  hs.OrthogonalPolynomial1D.jacobi.base = ...
    fullfile(hs.OrthogonalPolynomial1D.base, 'jacobi');
