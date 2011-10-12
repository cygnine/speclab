function[w] = constant_w(self,lambda,mu)
% constant_w -- The constant w for generalized Gegenbauer polynomials.
%
% w = constant_w(self, [[lambda=self.lambda, mu=self.mu]])
%
%     Returns the value of the constant w_{lambda, mu}, which is the inverse of
%     the integral of the function
%
%            f(x) = (1-x^2)^(lambda-1/2) * |x|^(2*mu),
%
%     over x \in [-1,1]. Therefore, w is the multiplicative normalization
%     constant to make f(x) unit-normalized in L^1 over [-1,1].

if nargin < 3
  mu = self.mu;
  if nargin < 2
    lambda = self.lambda;
  end
end

w = gammaln(lambda+mu+1)./gammaln(mu+1/2)./gammaln(lambda+1/2);
w = exp(w);
