function[a,b] = recurrence(self, n)
% recurrence -- Returns the recurrence coefficients for the polynomial family
%
% [alpha,beta] = recurrence(self, n)
%
%     Computes and returns the three-term recurrence coefficients alpha_m and
%     beta_m for each m in the array n. The outputs alpha and beta are of the
%     same size as n, and the input n array should contain non-negative
%     integers.
%
%     This function is the preferred method for constructing the recurrence
%     coefficients, and it calls self.standard_recurrence. The main utility of
%     this function is that it builds in any factors for
%     self.weight_normalization into the coefficient beta_0.

nsize = size(n);
n = n(:);

[a,b] = self.standard_recurrence(n);

% Build in weight normalization + affine mapping factor
K = self.scale_weight(1);
%b(n==0) = 1/sqrt(K)*b(n==0);
%b(n==0) = 1/K*b(n==0);
b(n==0) = K*b(n==0);

a = reshape(a, nsize);
b = reshape(b, nsize);
