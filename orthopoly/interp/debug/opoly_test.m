% Script to see how the opoly interpolation behaves

clear
close all
from speclab.orthopoly1d.interp import monomial_lu opoly_li_lu 
from speclab.orthopoly1d.interp import monomial_coeffs opoly_li_coeffs 
from labtools import typelatex as ltex

from speclab.orthopoly1d.jacobi import gauss_quadrature as gq
from speclab.orthopoly1d.jacobi import jacobi_poly
from speclab.monomials import multimonomial

% 2D
N = 28;
theta = [gq(N, 'alpha', 0, 'beta', 0), zeros([N 1])];
theta = gq(N, 'alpha', -1/2, 'beta', -1/2)*[1 0.6];
theta = rand([N 2])*2-1;

[l,W,p,k,u] = monomial_lu(theta);
[lo, Wo, po, ko, uo] = opoly_li_lu(theta);

f = zeros([N 1]);
tlagrange = ceil(N/2);
f(tlagrange) = 1;

coeffs = monomial_coeffs(l,W,p,k,u,2,f);
coeffso = opoly_li_coeffs(lo,Wo,po,ko,uo,2,f);

thetaf = [linspace(-1,1,101).' zeros([101 1])];
thetaf = linspace(-1,1,101);
[xf,yf] = meshgrid(thetaf, thetaf);
thetaf = [xf(:) yf(:)];
mf = multimonomial(thetaf, 0:size(coeffs,1)-1, 'dim', 2)*coeffs;
of = jacobi_poly(thetaf, 0:size(coeffso,1)-1, 'alpha', 0, 'beta', 0, 'dim', 2)*coeffso;

tshape = sqrt(size(thetaf,1));
tshape = [tshape tshape];
contour(reshape(thetaf(:,1), tshape), reshape(thetaf(:,2), tshape), ...
      reshape(mf, tshape), linspace(-0.5,1,30));
hold on; plot(theta(:,1), theta(:,2), 'k.', 'markersize', 10);
plot(theta(tlagrange,1), theta(tlagrange,2), 'kx', 'markersize', 15);
axis equal;
axis([-1 1 -1 1]);
set(gca, 'fontsize', 16, 'fontweight', 'b');
ltex(xlabel('$x_1$', 'fontsize', 16));
ltex(ylabel('$x_2$', 'fontsize', 16));
colormap gray

figure
contour(reshape(thetaf(:,1), tshape), reshape(thetaf(:,2), tshape), ...
      reshape(of, tshape), linspace(-0.5,1,30));
hold on; plot(theta(:,1), theta(:,2), 'k.', 'markersize', 10);
plot(theta(tlagrange,1), theta(tlagrange,2), 'kx', 'markersize', 15);
axis equal;
axis([-1 1 -1 1]);
ltex(xlabel('$x_1$', 'fontsize', 16));
ltex(ylabel('$x_2$', 'fontsize', 16));
set(gca, 'fontsize', 16, 'fontweight', 'b');
colormap gray
