function[hs] = handles__()
% [HS] = HANDLES__()
%
%     Returns directory pointers for common module

% This is by default
hs.base = fileparts(mfilename('fullpath'));

hs.debug.base = fullfile(hs.base,'debug');

hs.common.base = fullfile(hs.base,'common');
hs.fourier.base = fullfile(hs.base,'fourier');
  hs.fourier.eval.base = fullfile(hs.fourier.base, 'eval');
  hs.fourier.quad.base = fullfile(hs.fourier.base, 'quad');
  hs.fourier.weights.base = fullfile(hs.fourier.base, 'weights');
  hs.fourier.maps.base = fullfile(hs.fourier.base, 'maps');
  hs.fourier.connection.base = fullfile(hs.fourier.base, 'connection');
  hs.fourier.fft.base = fullfile(hs.fourier.base, 'fft');
hs.wiener.base = fullfile(hs.base,'wiener');
  hs.wiener.eval.base = fullfile(hs.wiener.base, 'eval');
  hs.wiener.weights.base = fullfile(hs.wiener.base, 'weights');
  hs.wiener.maps.base = fullfile(hs.wiener.base, 'maps');
hs.monomials.base = fullfile(hs.base,'monomials');
hs.NewtonPolynomials.base = fullfile(hs.base,'NewtonPolynomials');

hs.OrthogonalPolynomial1D.base = fullfile(hs.base,'OrthogonalPolynomial1D');
  hs.OrthogonalPolynomial1D.jacobi.base = ...
    fullfile(hs.OrthogonalPolynomial1D.base, 'jacobi');
  hs.OrthogonalPolynomial1D.jacobi.quad.base = ...
    fullfile(hs.OrthogonalPolynomial1D.jacobi.base, 'quad');
  hs.OrthogonalPolynomial1D.jacobi.eval.base = ...
    fullfile(hs.OrthogonalPolynomial1D.jacobi.base, 'eval');
  hs.OrthogonalPolynomial1D.jacobi.coefficients.base = ...
    fullfile(hs.OrthogonalPolynomial1D.jacobi.base, 'coefficients');
  hs.OrthogonalPolynomial1D.jacobi.connection.base = ...
    fullfile(hs.OrthogonalPolynomial1D.jacobi.base, 'connection');
