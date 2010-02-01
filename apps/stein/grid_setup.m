function[grid] = grid_setup(alpha,varargin)
% grid_setup -- Precomputes constants for online evaluation of Stein integrals
%
% grid = grid_setup(alpha,{interval=[0, 2*pi], N_s=10, N_theta=100})
%
%     Given a value of alpha, this function precomputes some quantities for the
%     second-difference Stein integral of order alpha. N_s and N_theta
%     represent the quadrature order in the s- and theta- directions,
%     respectively. The integral in the s-direction is truncated after a certain
%     number of loops based on the value of the Stein weight.

persistent fourier_gq jacobi_gq replicate_local_nodes ...
           strict_inputs
if isempty(strict_inputs)
  from labtools import strict_inputs
  from speclab.orthopoly1d.jacobi.quad import gauss_quadrature as jacobi_gq
  from speclab.fourier.quad import gauss_quadrature as fourier_gq
  from piecewise_interpolation.grid_tools import replicate_local_nodes
end

opt = strict_inputs({'interval', 'N_theta', 'N_s'}, {[0, 2*pi], 100, 10}, [], varargin{:});

temp.scale = diff(opt.interval)/(2*pi);
temp.shift = pi + opt.interval(1);

% First generate theta grid:
[grid.theta,grid.theta_w] = fourier_gq(opt.N_theta, temp);
dtheta = abs(grid.theta(2) - grid.theta(1));

temp.shift = 0; temp.scale = 1; temp.alpha = 0; temp.beta = 0;
% Now generate the s-grid that will be used on each subinterval between theta
% nodes:
[s,ws] = jacobi_gq(opt.N_s, temp);

% Now assume that theta is equispaced. Create a grid of s nodes that are
% spaced between each set of two theta nodes.
grid.global_s_nodes = replicate_local_nodes(s, [grid.theta; grid.theta(1) + 2*pi]);
grid.global_s_nodes(:,end) = mod(grid.global_s_nodes(:,end), 2*pi);
% (the last column straddles the interface)

grid.ws = ws/2*dtheta;

tol = 1e-10;
max_rotations = 1e5;
n_rotations = min(ceil(1/(2*pi*tol^(1/(2*alpha+1)))), max_rotations);

weights = repmat(grid.ws, [1, opt.N_theta]);
temp2 = mod(grid.global_s_nodes - grid.theta(1), 2*pi);
temp = 0*grid.global_s_nodes;
for q = 1:n_rotations
  temp = temp + 1./(temp2 + (q-1)*2*pi).^(2*alpha+1);
end
grid.weights = weights.*temp;

exponent = (2*alpha + 1);
if abs(exponent-2)<1e-10
  grid.b = pi/4;
elseif abs(exponent-3)<1e-10
  grid.b = log(2);
elseif abs(exponent-4)<1e-10
  grid.b = pi/3;
else
  % All hail Wolfram Alpha
  grid.b = 2^(exponent-5)*(2^exponent-8)*sin((pi*exponent)/2)*gamma(1-exponent);
end
grid.b = 1./(4*grid.b);

% WTF is this factor about?
grid.b = 2^(2*alpha-2)*grid.b;
