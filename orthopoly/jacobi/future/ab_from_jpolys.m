function[alphas,betas] = ab_from_jpolys(a,b,w,N,a0,b0);

% [alphas,betas] = ab_from_jpolys(a,b,w,N);
% Derives the N recurrence coefficients alphas and betas for orthogonal
% polynomials under a weight function over [-1,1] that is a linear combination 
% of Jacobi weight functions. 
% 
% Inputs a,b, and w all have length N. The total weight function is \sum_j w(j)
% weight^{a(j),b(j)}, where weight^{a,b} is the Jacobi weight function.
% 
% a0 and b0 are scalars that define the Jacobi polynomials which one starts from
% to derive the alphas and betas. (See chebyshev_modified_moments.m)
% 
% 20080522: acn

% Find moments nus
nus = zeros([1 2*N]);
M = length(a);

for q = 1:M;
  nus = nus + w(q)*calc_jacobi_moments(a(q),b(q),2*N,a0,b0);
end

% Get `known' recurrence coefficients
[as,bs] = jacobi_recurrence(2*N,a0,b0);

% Magic:
[alphas, betas] = chebyshev_modified_moments(nus,as,bs);
