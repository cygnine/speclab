function[x] = linspace(self, N)
% linspace -- Returns a vector of equispaced points on the interval
%
% x = linspace(self, N)
%
%     Returns a size-N vector of equispaced points on the interval (including
%     the endpoints).

x = linspace(self.interval(1), self.interval(2), N);
