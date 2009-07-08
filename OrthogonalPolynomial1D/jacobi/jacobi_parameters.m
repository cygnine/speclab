% Script that defaults the values Jacobi polynomial parameters if they
% aren't already defined.
% Interval is [shift-scale, shift+scale]
%
% 20080623 -- acn

% Defaults: alpha = 0, beta = 0

if ~exist('alpha','var')
  alpha = -1/2;
end
if ~exist('beta','var')
  beta = -1/2;
end
if ~exist('shift','var')
  shift = 0;
end
if ~exist('scale','var')
  scale = 1;
end
if ~exist('a','var');
  a = -scale + shift;
end
