function[h] = l2_norm(self,n)
% l2_norm -- Weighted l2 norm of polynomials
%
% h = l2_norm(self,n)
%
%     Returns the weighted L^2 norm of the polynomials of degree n. This is a
%     function both of the weight normalization and function normalization of
%     the class instance.

% First get weight scaling:
K = self.map_to_domain.A*self.scale_weight(1);

% Now get function scalings:
h = self.scale_functions(ones(size(n)), n);

h = h.*sqrt(K);
