function[grid] = gauss_grid_setup(alpha, varargin)
% gauss_grid_setup -- Precomputes constants for a Gauss grid 
%
% grid = grid_setup(alpha, {theta=linspace(0, 2*pi, 100),  N_s=100, N_p = 3})
%
%     Given a value of alpha, this function precomputes some quantities for the
%     second-difference Stein integral of order alpha. N_s 
%     represents the quadrature order in the s-direction.
%     The integral in the s-direction uses semi-infinite quadrature with
%     the Moebius mapping taking [-1,1] ---> [0, \infty). The input N_p
%     denotes the order of approximation on each theta-cell. (I.e. it's a
%     finite-element representation.) 
%
%     Theta represents the cell vertices of the finite element grid. An order
%     N_p element is formed for each cell specified by theta. 

persistent fourier_gq jacobi_gq replicate_local_nodes ...
           strict_inputs stein_proportionality spdiag
if isempty(strict_inputs)
  from labtools import strict_inputs spdiag
  from speclab.orthopoly1d.jacobi.quad import gauss_quadrature as jacobi_gq
  from speclab.fourier.quad import gauss_quadrature as fourier_gq
  from piecewise_interpolation.grid_tools import replicate_local_nodes
  from speclab.apps.stein import second_difference_sobolev_proportionality as stein_proportionality
end

opt = strict_inputs({'theta', 'N_s', 'N_p'}, ...
                    {linspace(0, 2*pi, 100), 100, 3}, [], varargin{:});

% Grid on each theta-element:
[grid.standard_element.r, grid.standard_element.w] = jacobi_gq(opt.N_p, 'alpha', 0, 'beta', 0);

% *Sigh* I haven't implemented this part of Speclab yet...time to do things
% the hard way:
[grid.s, grid.w_s] = jacobi_gq(opt.N_s, 'alpha', 0, 'beta',0);
grid.s = (1 + grid.s)./(1-grid.s);
grid.w_s = 1/2*grid.w_s.*(grid.s+1).^2./grid.s.^(2*alpha+1);
grid.N_s = opt.N_s;
% Ok, now (grid.s, grid.w_s) is a quadrature rule integrating against weight
% 1/s^(2*alpha+1). (Note the extreme suckage if alpha < 1/2.)

% One more step: scale it so that one half rotation is encompassed in the
% Jacobi [-1, 0] interval quadrature. 
%scale = 0.5;
%grid.s = grid.s*scale;
%grid.w_s = grid.w_s/scale^(2*alpha);

% Use theta, N_p to generate a finite element mesh
grid.theta_grid = replicate_local_nodes(grid.standard_element.r, opt.theta);
grid.w_theta = repmat(grid.standard_element.w, [1 size(grid.theta_grid,2)]);
grid.theta_scales = diff(opt.theta(:))/2;
grid.w_theta = grid.w_theta*spdiag(grid.theta_scales);
grid.N_p = opt.N_p;

grid.b = stein_proportionality(alpha);
