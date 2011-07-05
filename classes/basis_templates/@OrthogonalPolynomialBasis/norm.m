function[h] = norm(self,n)
% norm -- (native) Hilbert norm of polynomials
%
% h = norm(self,n)
%
%     Returns the weighted L^2 norm of the polynomials of index n. This is a
%     function both of the weight normalization and function normalization of
%     the class instance.

% First get weight scaling:
K = self.map_to_domain.A*self.scale_weight(1);

% Now get function scalings:
h = self.scale_functions(ones(size(n)), n);

h = h.*sqrt(K);
