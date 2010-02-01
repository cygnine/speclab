% Script to test accuracy of Stein integrals

clear 
close all

from speclab.apps import stein
from speclab.fourier.quad import gauss_quadrature as fourier_gq
from speclab.fourier.eval import fseries dfseries
from speclab.common import integer_range

fs = {@(x) sin(x), @(x) cos(x), @(x) sin(2*x), @(x) cos(2*x), ...
      @(x) sin(3*x), @(x) sin(4*x), @(x) sin(301*x)};

qs = 1:500;

Is = zeros([length(fs) 1]);
alpha = 1.5;

Js = zeros([length(fs) 1]);
N_theta = 25;
N_s = 20;
Nq = N_theta*N_s + N_theta;

% Use generic Fourier Series for comparison
temp.shift = pi; temp.scale = 1;
[theta,w] = fourier_gq(Nq, temp);
vandermonde = fseries(theta, integer_range(Nq), temp);
ks = integer_range(Nq);
dmatrix = vandermonde*spdiags((i*ks).^(alpha), 0, Nq, Nq)*inv(vandermonde);
%dvandermonde = dfseries(theta, integer_range(Nq),temp)*inv(vandermonde);
grid = stein.grid_setup(alpha, 'N_theta', 50, 'N_s', 10);

for q = qs
  f = @(x) sin(q*x);
  Is(q) = stein.stein_integral_d2_driver(f, grid);
  %Js(q) = sum(w.*abs((dvandermonde*f(theta)).^2));
  Js(q) = sum(w.*abs((dmatrix*f(theta)).^2));
end

exact_Is = pi*qs(:).^(2*alpha);

plot(qs, Is, 'r.', qs, Js, 'b.', qs, exact_Is, 'k');
