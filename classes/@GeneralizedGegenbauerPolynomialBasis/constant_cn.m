function[c] = constant_cn(self,n,lambda,mu)
% constant_cn -- The constant c_n for generalized Gegenbauer polynomials
%
% c = constant_cn(self, n, [[lambda=self.lambda, mu=self.mu]])
%
%     Returns an array c of the same size of n that evaluates the generalized
%     Gegenbauer constants c_n(lambda,mu). If lambda and mu are given as
%     arrays, they must be of the same size as n.
%
%     If G(.) is the Gamma function, the constant is given by 
%
%                           G(mu+1/2) * G(lambda+1/2) * G(n+lambda+mu) * G(n+1) * (2*n+lambda+mu)
%     (c_n(lambda,mu))^2 =  ---------------------------------------------------------------------
%                           G(lambda+mu+1) * G(n+mu+1/2) * G(n+lambda+1/2)
%
%     This si the normalization constant connecting the classical Jacobi
%     polynomials to the orthonormal generalized Gegenbauer polynomials.

Nsize = size(n);
if nargin < 4
  mu = ones(Nsize)*self.mu;
  if nargin < 3
    lambda = ones(Nsize)*self.lambda;
  end
end

c = ones(Nsize);

%c = gammaln(mu+1/2).*gammaln(lambda+1/2)./gammaln(mu+lambda+1/2);
c = c.*gammaln(n+lambda+mu).*gammaln(n+1)./gammaln(n+mu+1/2)./...
       gammaln(n+lambda+1/2);

c = exp(c/2);
c = c.*sqrt(2*n+lambda+mu);
