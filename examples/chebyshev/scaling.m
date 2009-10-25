%%% Example script for Chebyshev approximation
%%% Scaling: how to change intervals of approximation

% In Matlab helpstrings for functions, optional arguments are given in curly
% brackets.  E.g., the helpstring declaration
%
% [output] = myfun(in1, in2, {opt1=value1, opt2=value2})
%
% means that myfun takes two required arguments and two optional arguments; the
% names of the optional variables are 'opt1' and 'opt2', and they default to
% value1 and value2, respectively.
%
% To give optional arguments into *any* routines in speclab, you should use the
% ('key', value) input syntax. E.g., to call myfun and specify opt2 to be 3, you
% should call
%
% >> myfun(whatever1, whatever2, 'opt2', 3)
%
% To specify both parameters, you can call
%
% >> myfun(whatever1, whatever2, 'opt1', something, 'opt2', otherthing)
%      ---- OR ----
% >> myfun(whatever1, whatever2, 'opt2', otherthing, 'opt1', something)
%
% I.e., the order of optional inputs doesn't matter. There is an alternative to
% specifying ('key', value) pairs; see below.
%
% Any optional arguments that aren't recognized are retained in each function,
% but are ignored. 

clear;
global packages;
speclab = packages.speclab;
cheb = speclab.orthopoly1d.jacobi;
common = packages.common;

%% Manually set affine scale and shift parameters.
scale = 3;
shift = -1; 
% This sets the interval of approximation to be [-1*scale+shift, 1*scale+shift]
% = [-4, 2]. Let's do Gauss-Lobatto quadrature to highlight this.
N = 50;
[r,w] = cheb.quad.gauss_lobatto_quadrature(N,'scale',scale,'shift',shift);

fprintf('The smallest nodal value is %1.2f and the largest nodal value is %1.2f\n', ...
  min(r), max(r));

%% Simplification 1: bundling optional arguments
% Specifying optional arguments in ('key', value) pairs is cumbersome for large
% numbers of optional arguments. To ameliorate this, you can bundle key/value
% pairs into a struct, in which case you can bypass the onerous ('key', value)
% input scheme and just input the whole struct. The struct must be the only
% optional argument given.
map.scale = 3;
map.shift = -1;
[r,w] = cheb.quad.gauss_lobatto_quadrature(N,map);
fprintf('The smallest nodal value is %1.2f and the largest nodal value is %1.2f\n', ...
  min(r), max(r));

% NOTE: you *cannot* mix and match struct+('key',value) input schemes. I.e., the
% following will *not* work:
% >> opt.opt1 = 3;
% >> myfun(blah1, blah2, 'opt2', something, opt);

%% Simplification 2: have speclab calculate affine parameters for you
% Most of the time, you want to specify the interval of approximation, not the
% affine parameters explicitly. Speclab can handle this for you.
map = cheb.affine_scaling([-4,2]);
[r,w] = cheb.quad.gauss_lobatto_quadrature(N,map);
fprintf('The smallest nodal value is %1.2f and the largest nodal value is %1.2f\n', ...
  min(r), max(r));

% Many files in speclab have scale and shift optional arguments. You can supply
% them to any function taking these inputs. (Usually evaluation and quadrature
% routines, and some operator/coefficient/matrix/fft routines). See the
% helpstring for each function to determine if scaling and shifting is part of
% the inputs.
