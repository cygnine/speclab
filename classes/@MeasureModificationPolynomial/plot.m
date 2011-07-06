function[varargout] = plot(self, varargin)
% plot -- Plots the MeasureModificationPolynomial over the defined interval
%
% plot(self, ...)

x = self.interval.linspace(100);
[varargout{1:nargout}] = plot(x, self.evaluate(x), varargin{:});
