function[a,b] = standard_recurrence(self, n)
% standard_recurrence -- Returns the 'standard' three-term recurrence coefficients
%
% [alpha,beta] = standard_recurrence(n)
%
%     Computes and returns the recurrence coefficients alpha_m and beta_m for
%     each m in the array n. 
%
%     On the standard interval self.standard_domain, all orthogonal polynomial
%     systems with orthonormal elements satisfy the recurrence relation:
%
%       sqrt(beta_{n+1}) p_{n+1} = (x - alpha_n) p_n - sqrt(beta_n) p_{n-1},
%       p_0 = 1/sqrt(beta_0),
%       p_{-1} = 0
%
%     The coefficients alpha and beta for each index in n are returned, and the
%     output arrays have the same size as the input array n.

nsize = size(n);
n = n(:);

if not(isempty(self.recurrence_handle))
  [a,b] = self.recurrence_handle(n);
else
  error('You have not defined the "recurrence_handle" property of this class instance');
end

a = reshape(a, nsize);
b = reshape(b, nsize);
