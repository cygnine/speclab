%%% Example script for Legendre approximation
%%% Approximation convergence: tests convergence of spectral approximations

% Test a hallmark of `spectral' methods: exponential convergence for analytic
% functions. 
clear
global packages;
exloglog = packages.common.exloglog;
typelatex = packages.common.typelatex;
leg = packages.speclab.orthopoly1d.jacobi;
opt.alpha = 0;
opt.beta = 0;

% Maximum quadrature size
Nq = 500;
% Maximum # of modes
N = 150;
[r,w] = leg.quad.gauss_quadrature(Nq,opt);

f = @(r) sin(45*r);
df = @(r) 45*cos(45*r);

% Compute all modes
fr = f(r);
dfr = df(r);
v = leg.eval.eval_jacobi_poly(r,0:(N-1),opt);
modes = v'*(w.*fr);
diffmodes = leg.operators.stiffness_operator(modes,opt);

% Initialize error vectors
Nplot = 10:(N-1);
f_l2_err = zeros([1,length(Nplot)]);
f_linf_err = f_l2_err;
df_l2_err = f_l2_err;
df_linf_err = f_l2_err;
qcount = 1;

% Do all the hard work
for N=Nplot
  temp = modes;
  temp((N+1):end) = 0;
  approx = v*temp;
  f_l2_err(qcount) = sqrt(sum(w.*(approx-fr).^2));
  f_linf_err(qcount) = max(abs(approx-fr)); 

  temp = diffmodes;
  temp((N+1):end) = 0;
  approx = v*temp;
  df_l2_err(qcount) = sqrt(sum(w.*(approx-dfr).^2));
  df_linf_err(qcount) = max(abs(approx-dfr));

  qcount = qcount + 1;
end

% The remainder is plotting and plot formatting
figure;
subplot(2,2,1)
exloglog(Nplot,f_l2_err);
typelatex(xlabel('$N$'));
typelatex(ylabel('$L^2$ error'));
temp = axis; 
axis([min(Nplot), max(Nplot), 0.1*min(f_l2_err), 10*max(f_l2_err)]);
temp = axis;
temp = [floor(log10(temp(3))), ceil(log10(temp(4)))];
labels = 10.^(temp(1):2:temp(2));
set(gca, 'YTick', labels);
set(gca, 'YTickLabel', labels);

subplot(2,2,2)
exloglog(Nplot,f_linf_err);
typelatex(xlabel('$N$'));
typelatex(ylabel('$L^\infty$ error'));
temp = axis; 
axis([min(Nplot), max(Nplot), 0.1*min(f_linf_err), 10*max(f_linf_err)]);
temp = axis;
temp = [floor(log10(temp(3))), ceil(log10(temp(4)))];
labels = 10.^(temp(1):2:temp(2));
set(gca, 'YTick', labels);
set(gca, 'YTickLabel', labels);

subplot(2,2,3)
exloglog(Nplot,df_l2_err);
typelatex(xlabel('$N$'));
typelatex(ylabel('$L^2$ derivative error'));
temp = axis; 
axis([min(Nplot), max(Nplot), 0.1*min(df_l2_err), 10*max(df_l2_err)]);
temp = axis;
temp = [floor(log10(temp(3))), ceil(log10(temp(4)))];
labels = 10.^(temp(1):2:temp(2));
set(gca, 'YTick', labels);
set(gca, 'YTickLabel', labels);

subplot(2,2,4)
exloglog(Nplot,df_linf_err);
typelatex(xlabel('$N$'));
typelatex(ylabel('$L^\infty$ derivative error'));
temp = axis; 
axis([min(Nplot), max(Nplot), 0.1*min(df_linf_err), 10*max(df_linf_err)]);
temp = axis;
temp = [floor(log10(temp(3))), ceil(log10(temp(4)))];
labels = 10.^(temp(1):2:temp(2));
set(gca, 'YTick', labels);
set(gca, 'YTickLabel', labels);
