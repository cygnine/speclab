%%% Example script for Legendre approximation
%%% Resolution scaling: doing some advanced affine scaling methods

% The 'affine_scaling' function in speclab.orthopoly.jacobi also packages
% `fractional resolution scaling'. The resolution each collection of basis
% functions can attain can be characterized using their 'canonical' nodal
% distributions. For orthogonal polynomials, this canonical set of modes is the
% Gauss quadrature nodal set. We can use this characterization to decide how to
% scale the basis if I want, say, 80% of the resolution inside some given
% interval. 

clear
import speclab.orthopoly1d.jacobi as leg
%leg = packages.speclab.orthopoly1d.jacobi;
opt.alpha = 0; 
opt.beta = 0;

% The following three calls to affine scaling do the same thing. Note that since
% we're taking Legendre polys as a special case of Jacobi polys, it is necessary
% to specify alpha=beta=0.
interval = [0,1];
  % Put 100% of the resolution inside the specified interval
  leg.affine_scaling(interval, 'resolution_fraction', 1, 'alpha', 0, 'beta', 0);
  % Put 100% of the resolution for 100 modes inside the specified interval
  leg.affine_scaling(interval, 'N', 100, 'resolution_fraction', 1, 'alpha', 0, 'beta', 0);
  % Put 80% of the resolution for an unspecified # of modes inside interval
  % This is an ill-posed request since N is unspecified
  leg.affine_scaling(interval, 'resolution_fraction', 0.8, 'alpha', 0, 'beta', 0);

% The following call says 'return the scale/shift parameters to put 73% of the
% resolution for 250 modes inside [0,pi]':
opt.resolution_fraction = 0.73;
opt.N = 250;
map = leg.affine_scaling([0,pi], opt);
opt.shift = map.shift; opt.scale = map.scale;

% You can check that the shift and scale returned are different for any given
% input value of N if you fix resolution_fraction.
% Check that what we asked for is what we got:
[r,w] = leg.quad.gauss_quadrature(opt.N,opt);

% Now we don't know what the interval of approximation is, we only know the
% affine shift and scale parameters. While it's a quick arithmetic calculation
% to determine the interval of approximation, just let a computer do it for you.
str = ['For 73 percent of a 250-point Legendre resolution to lie inside '...
       '[0,pi],\nthe interval of approximation must be [%1.3f, %1.3f]\n'];
str2 = ['With this interval of approximation, a ratio %1.4f of the 250\nGauss nodes lie' ...
        ' inside [0,pi].\n'];
approx_interval = leg.interval(opt);
fprintf(str,approx_interval(1),approx_interval(2));
fprintf(str2,sum(r>=0 & r<=pi)/opt.N);
