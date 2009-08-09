function[hs,pathadditions] = handles__()
% [HS,PATHADDITIONS] = HANDLES__()
%
%     Returns directory pointers for common module in HS. PATHADDITIONS is a
%     cell array with a string in each element indicated paths to add to the
%     global path structure. 

% This is by default
hs.base = fileparts(mfilename('fullpath'));

hs.debug.base = fullfile(hs.base,'debug');

hs.common.base = fullfile(hs.base,'common');

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
  hs.OrthogonalPolynomial1D.jacobi.weights.base = ...
    fullfile(hs.OrthogonalPolynomial1D.jacobi.base, 'weights');
  hs.OrthogonalPolynomial1D.jacobi.operators.base = ...
    fullfile(hs.OrthogonalPolynomial1D.jacobi.base, 'operators');

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

pathadditions = cell(0);
pathadditions{end+1} = fullfile(hs.base, 'classes');
