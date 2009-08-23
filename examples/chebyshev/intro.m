%%% Example script for Chebyshev approximation
%%% Introduction: calling syntax of basic routines

% In most of the Chebyshev (Jacobi) routines, the optional parameters 'scale'
% and 'shift' are 1 and 0 by default, respectively. If you change these
% parameters, the approximation will be valid on the interval [-scale+shift,
% scale+shift]. This corresponds to a mapping r -----> s = scale*r + shift.

clear
% The following two lines extract the 'module' of functions we'll need
global handles;
cheb = handles.speclab.orthopoly1d.jacobi;
common = handles.common;
% We've called 'cheb' the jacobi module. Unless you say otherwise, speclab
% assumes the Jacobi polynomial you want to evaluate is the Chebyshev kind. To
% see how to specify a different Jacobi polynomial type, see the jacobi examples
% folder. 

%% Construct a chebyshev Vandermonde matrix on the default [-1,1] interval.
N = 50;  % or whatever
[r,w] = cheb.quad.gauss_quadrature(N);   % N-point chebyshev-gauss quadrature
vandermonde = cheb.eval.eval_jacobi_poly(r,0:(N-1)); % The L2-normalized polys

% To get help on any function, use the function handles.helper. You need to know
% the full tree structure of the file you want to know information about. For
% example, print the helpstring for cheb.quad.gauss_quadrature:
fprintf('Helpstring for Jacobi Gauss quadrature routine:\n');
handles.helper('speclab','orthopoly1d','jacobi','quad','gauss_quadrature')
% Alternatively, you can also cd into the appropriate directory and use Matlab's
% help function, but this latter approach requires more typing.

% Define a function
f = @(r) sin(5*(r-3).^2);

% Determine an interpolant by inverting the Vandermonde matrix
modal_coefficients_1 = inv(vandermonde)*f(r);

% Better way: use the quadrature rule. The mass matrix is the identity since by
% default the vandermonde matrix is the evaluation of L^2 normalized polynomials
modal_coefficients_2 = vandermonde'*(f(r).*w);

% The following is effectively the L2 error
err = norm(modal_coefficients_1 - modal_coefficients_2);
fprintf('Error between quadrature rule and matrix inversion is %1.4e\n', err);

%% Plotting interpolants
% Let's see how accurate the interpolant is with a plot
r_refined = linspace(-1,1,1e3).'; % *must* be a column vector...otherwise error
vandermonde_refined = cheb.eval.eval_jacobi_poly(r_refined,0:(N-1));

common.exsemilogy(r_refined,abs(f(r_refined) - vandermonde_refined*modal_coefficients_2));
common.typelatex(xlabel('$\mathbf{r}$'));
common.typelatex(ylabel('$|f(r) -- f_N(r)|$'));

% If you *really* wanted to, you could use the monic polys instead:
vandermonde = cheb.eval.eval_jacobi_poly(r,0:(N-1), 'normalization', 'monic');
modal_coefficients_1 = inv(vandermonde)*f(r);

% now the mass matrix isn't the identity anymore. Figure it out explicitly:
mass_weights = diag(vandermonde'*spdiags(w,0,N,N)*vandermonde);
modal_coefficients_2 = (vandermonde'*(f(r).*w))./mass_weights;

err = norm((modal_coefficients_1 - modal_coefficients_2).*sqrt(mass_weights));
fprintf('Error between quadrature rule and matrix inversion (monic polys) is %1.4e\n', err);

% NOTE: check cond(vandermonde) for monic polys...it's horrible. Use quadrature
% instead.
