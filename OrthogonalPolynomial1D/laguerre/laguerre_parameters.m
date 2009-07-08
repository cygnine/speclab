% Script that defaults the values Laguerre polynomial parameters if they
% aren't already defined.
%
% 20080623 -- acn

% Defaults: alpha = 0, shift = 0;

if ~exist('alpha','var')  % Exponent on weight x
  alpha = 0;
end
if ~exist('shift','var')  % Rigid shift of domain
  shift = 0;
end
if ~exist('scale','var')  % Linear scaling of domain
  scale = 1;
end
if ~exist('xa','var'); % point for Gauss-Radau quadrature
  xa = 0;
end
