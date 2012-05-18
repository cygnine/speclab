function[w] = weight(self, x)
% weight -- Evaluates weight function on the unit disc
%
% w = weight(self, x)
%
%     Given 2D coordinates x (a two-column vector) this function evaluates the
%     weight function on the unit ball evaluated at these points. The weight
%     function (modulo user-specified normalization) is
%
%         w(x,y) = (1 - x^2 - y^2)^(self.mu-1/2)

assert(size(x,2)==2, 'Input x must have two columns');

r = self.map_to_standard_domain(x);
w = (1 - r(:,1)^.2 - r(:,2).^2).^(self.mu-1/2);
w = self.scale_weight(w);
