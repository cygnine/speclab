%%% Example script for Wiener approximation
%%% Plotting functions: plots Wiener basis functions

clear
global packages;
wiener = packages.speclab.wiener;
explot = packages.labtools.explot;
ltex = packages.labtools.typelatex;

opt.s = 3;
N = 500;
[x,w] = wiener.quad.pi_gauss_quadrature(N,opt);
% Some notes:
% - with no optional inputs, s=1, which is the regular Fourier case.
% - you can still scale/shift with the Wiener basis, using the function
%   affine_scaling, but scale/shift optional input do not change the interval of
%   approximation: it is always the entire real line. However, scale/shift
%   parameters do change the resolution of the basis functions, and it's
%   useful for example to say "I want all the resolution for an s=5 N=43 Wiener
%   expansion to lie inside the interval [-3,4]".
% - "pi_gauss_quadrature" : since the Wiener functions are weighted versions of
%   the Fourer series, we don't use the canonical gauss_quadrature, but instead
%   a weighted version of it, called pi_gauss_quadrature. If you want to use the
%   *unweighted* Wiener functions, then you can use the regular
%   gauss_quadrature.
% - Due to the definition of how the Wiener functions are weighted, then if s is
%   an integer, the Gauss quadrature rule will actually successfully integrate
%   the Wiener functions. However, some over-integration is necessary.

x_plot = linspace(-3,3,1000).';
ws = wiener.eval.wiener_function(x_plot,0:4,opt);

subplot(2,1,1);
explot(x_plot, real(ws)); 
ltex(xlabel('$x$')); 
ltex(ylabel(['$Re\{\phi_k^{(' num2str(opt.s) ')}\}$']));
temp = axis; axis([min(x_plot), max(x_plot), temp(3:4)]);
ltex(title('Wiener function plots, $k=0,1,2,3,4$'));

subplot(2,1,2);
explot(x_plot, imag(ws)); 
ltex(xlabel('$x$')); 
ltex(ylabel(['$Im\{\phi_k^{(' num2str(opt.s) ')}\}$']));
temp = axis; axis([min(x_plot), max(x_plot), temp(3:4)]);

