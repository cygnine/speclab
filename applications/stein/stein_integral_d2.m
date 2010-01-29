function[I] = stein_integral_d2(f, alpha)
% stein_integral_d2(f, alpha)
%
% I = stein_integral_d2(f, alpha)
%
%     Given a function f : [0, 2*pi] -----> R that is periodic, this function
%     computes:
%
%                   2\pi    \infty
%                 /        / 
%                 |        |  | f(theta-s) + f(theta+s) - 2*f(theta) |^2
%     b(alpha) *  |        | -------------------------------------------  d(s) d(theta)
%                 |        |               s^{2*alpha+1}
%                 /        / 
%                 0        -\infty 
%
%     It uses Fourier quadrature for the theta integral and Gauss-Legendre
%     quadrature for the s integral. b(alpha) is the factor that makes the
%     scales the integral to equal the L^2 Sobolev(alpha) seminorm.

persistent grid_setup driver
if isempty(grid_setup)
  from speclab.applications.stein import grid_setup
  from speclab.applications.stein import stein_integral_d2_driver as driver
end

grid = grid_setup(alpha, 'N_theta', 50, 'N_s', 10);
I = driver(f, grid);
