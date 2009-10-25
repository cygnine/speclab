function[hs,pathadditions,name] = packages__(varargin)
% packages__ -- constructs path tree for the speclab module
%
% [hs,pathadditions,name] = packages__()
%
%     Returns directory pointers for common module in hs. pathadditions is a
%     cell array with a string in each element indicated paths to add to the
%     global path structure. name is the name of this module.

name = 'speclab';

% This is by default
hs.base = fileparts(mfilename('fullpath'));

hs.common.base = fullfile(hs.base,'common');

hs.orthopoly1d.base = fullfile(hs.base,'orthopoly1d');
  hs.orthopoly1d.jacobi.base = fullfile(hs.orthopoly1d.base, 'jacobi');
    hs.orthopoly1d.jacobi.quad.base = fullfile(hs.orthopoly1d.jacobi.base, 'quad');
    hs.orthopoly1d.jacobi.eval.base = fullfile(hs.orthopoly1d.jacobi.base, 'eval');
    hs.orthopoly1d.jacobi.coefficients.base = ...
      fullfile(hs.orthopoly1d.jacobi.base, 'coefficients');
    hs.orthopoly1d.jacobi.connection.base = ...
      fullfile(hs.orthopoly1d.jacobi.base, 'connection');
    hs.orthopoly1d.jacobi.weights.base = ...
      fullfile(hs.orthopoly1d.jacobi.base, 'weights');
    hs.orthopoly1d.jacobi.operators.base = ...
      fullfile(hs.orthopoly1d.jacobi.base, 'operators');
    hs.orthopoly1d.jacobi.fft.base = fullfile(hs.orthopoly1d.jacobi.base, 'fft');

  hs.orthopoly1d.hermite.base = fullfile(hs.orthopoly1d.base, 'hermite');
    hs.orthopoly1d.hermite.quad.base = fullfile(hs.orthopoly1d.hermite.base, 'quad');
    hs.orthopoly1d.hermite.eval.base = fullfile(hs.orthopoly1d.hermite.base, 'eval');
    hs.orthopoly1d.hermite.coefficients.base = ...
      fullfile(hs.orthopoly1d.hermite.base, 'coefficients');
    hs.orthopoly1d.hermite.weights.base = ...
      fullfile(hs.orthopoly1d.hermite.base, 'weights');

  hs.orthopoly1d.laguerre.base = fullfile(hs.orthopoly1d.base, 'laguerre');
    hs.orthopoly1d.laguerre.quad.base = fullfile(hs.orthopoly1d.laguerre.base, 'quad');
    hs.orthopoly1d.laguerre.eval.base = fullfile(hs.orthopoly1d.laguerre.base, 'eval');
    hs.orthopoly1d.laguerre.coefficients.base = ...
      fullfile(hs.orthopoly1d.laguerre.base, 'coefficients');
    hs.orthopoly1d.laguerre.weights.base = ...
      fullfile(hs.orthopoly1d.laguerre.base, 'weights');
    hs.orthopoly1d.laguerre.connection.base = ...
      fullfile(hs.orthopoly1d.laguerre.base, 'connection');
    hs.orthopoly1d.laguerre.operators.base = ...
      fullfile(hs.orthopoly1d.laguerre.base, 'operators');

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

hs.newton_polynomials.base = fullfile(hs.base,'newton_polynomials');

hs.debug.base = fullfile(hs.base,'debug');
hs.examples.base = fullfile(hs.base,'examples');
  hs.examples.chebyshev.base = fullfile(hs.examples.base, 'chebyshev');
  hs.examples.legendre.base = fullfile(hs.examples.base, 'legendre');
  hs.examples.jacobi.base = fullfile(hs.examples.base, 'jacobi');
  hs.examples.fourier.base = fullfile(hs.examples.base, 'fourier');
  hs.examples.wiener.base = fullfile(hs.examples.base, 'wiener');

pathadditions = cell(0);
pathadditions{end+1} = fullfile(hs.base, 'classes');
