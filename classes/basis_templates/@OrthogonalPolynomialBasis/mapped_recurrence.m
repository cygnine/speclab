function[alpha,beta,gamma] = mapped_recurrence(self, n)
% mapped_recurrence -- Recurrence coefficients for mapped functions
%
% [a,b,c] = mapped_recurrence(self, n)
%
%     For the particular normalization of the instance, returns the vector of
%     coefficients a,b,c such that 
%
%      a_{n+1} p_{n+1} = (x + b_{n+1}) p_n + c_{n+1} p_{n-1}.
%      a_0 p_0 = 1
%      p_{-1} = 0
%
%     Note that the coefficients in this relation are different from those in
%     the orthonormal and/or monic recurrence relation.

N = max(n(:));
ns = (0:N).';

% Old recurrence + leading coefficients (compared to orthonormal normalization)
% needed
[a,b] = self.recurrence(ns);
c = self.scale_functions(ones([1 N+1]), ns).';

alpha = zeros([N+1 1]);
beta = zeros([N+1 1]);
gamma = zeros([N+1 1]);

% New recurrence in terms of old coefficients on standard domain
alpha(2:end) = sqrt(b(2:end)).*c(1:end-1)./c(2:end);
beta(2:end) = -a(1:end-1);
gamma(3:end) = -sqrt(b(2:end-1)).*c(2:end-1)./c(1:end-2);
alpha(1) = c(1)/sqrt(b(1));

% Pluck out appropriate values of n
alpha = alpha(n+1);
beta = beta(n+1);
gamma = gamma(n+1);

% Assuming functions aren't scaled after mapped, this maps the recurrence
% relation to the physical domain
alpha = alpha./self.map_to_standard_domain.A;
beta = (self.map_to_standard_domain.b + beta)./self.map_to_standard_domain.A;
gamma = gamma./self.map_to_standard_domain.A;
