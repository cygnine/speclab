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

hs.orthopoly1d.base = fullfile(hs.base,'orthopoly1d');
  hs.orthopoly1d.jacobi.base = ...
    fullfile(hs.orthopoly1d.base, 'jacobi');
  hs.orthopoly1d.jacobi.quad.base = ...
    fullfile(hs.orthopoly1d.jacobi.base, 'quad');
  hs.orthopoly1d.jacobi.eval.base = ...
    fullfile(hs.orthopoly1d.jacobi.base, 'eval');
  hs.orthopoly1d.jacobi.coefficients.base = ...
    fullfile(hs.orthopoly1d.jacobi.base, 'coefficients');
  hs.orthopoly1d.jacobi.connection.base = ...
    fullfile(hs.orthopoly1d.jacobi.base, 'connection');
  hs.orthopoly1d.jacobi.weights.base = ...
    fullfile(hs.orthopoly1d.jacobi.base, 'weights');
  hs.orthopoly1d.jacobi.operators.base = ...
    fullfile(hs.orthopoly1d.jacobi.base, 'operators');
  hs.orthopoly1d.jacobi.jfft.base = ...
    fullfile(hs.orthopoly1d.jacobi.base, 'jfft');

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
  hs.wiener.quad.base = fullfile(hs.wiener.base, 'quad');
  hs.wiener.maps.base = fullfile(hs.wiener.base, 'maps');
  hs.wiener.fft.base = fullfile(hs.wiener.base, 'fft');
  hs.wiener.operators.base = fullfile(hs.wiener.base, 'operators');
  hs.wiener.matrices.base = fullfile(hs.wiener.base, 'matrices');
  hs.wiener.coefficients.base = fullfile(hs.wiener.base, 'coefficients');

hs.monomials.base = fullfile(hs.base,'monomials');

hs.NewtonPolynomials.base = fullfile(hs.base,'NewtonPolynomials');

pathadditions = cell(0);
pathadditions{end+1} = fullfile(hs.base, 'classes');
