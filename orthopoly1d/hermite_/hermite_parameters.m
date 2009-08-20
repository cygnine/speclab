% Script that defaults the values Hermite polynomial/function parameters if they
% aren't already defined.
%
% 20080623 -- acn

% Defaults: mu = 0, shift = 0, scale = 1

if ~exist('mu', 'var')
  mu = 0;
end
if ~exist('shift','var')
  shift = 0;
end
if ~exist('scale','var')
  scale = 1;
end
