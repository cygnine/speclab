function[A,B,C] = mapped_recurrence(self, n)
% mapped_recurrence -- Recurrence coefficients for mapped functions
%
% [A,B,C] = mapped_recurrence(self, n)
%
%     For the particular normalization of the instance, returns the vector of
%     coefficients a,b,c such that 
%
%      p_0 = B_0, A_0 = C_0 = 0
%      x p_1 = A_1 p_1 + B_1 p_0,  C_0 = 0
%      x p_{n+1} = A_{n+1} p_{n+1} + B_{n+1} p_n + C_{n+1} p_{n-1},   (n > 0)
%
%     This fucntion returns the coefficients {A_n, B_n, C_n} in arrays for each
%     integer in the input array n.

%      a_{n+1} p_{n+1} = (x + b_{n+1}) p_n + c_{n+1} p_{n-1}.
%      a_0 p_0 = 1
%      p_{-1} = 0

nsize = size(n);
n = n(:);
N = max(n);
ns = (0:N).';

% Old recurrence + leading coefficients (ratio to orthonormal normalization)
% needed
[a,b] = self.recurrence(ns);
kappa = self.scale_functions(ones([1 N+1]), ns).';

A = zeros([N+1 1]);
B = zeros([N+1 1]);
C = zeros([N+1 1]);

% Basic recurrence for orthonormal normalization, standard interval:
B(1) = 1/sqrt(b(1));
A(2:end) = sqrt(b(2:end));
B(2:end) = a(1:end-1);
C(3:end) = sqrt(b(2:end-1));

% First deal with any re-normalizations
% For p_0:
B(1) = kappa(1)*B(1);
% For the rest:
A(2:end) = (kappa(1:end-1)./kappa(2:end)).*A(2:end);
C(3:end) = (kappa(2:end-1)./kappa(1:end-2)).*C(3:end);

% Now deal with any mappings: 
% x = alpha*r + beta
alpha = self.map_to_domain.A;
beta = self.map_to_domain.b;

A(2:end) = alpha*A(2:end);
C(3:end) = alpha*C(3:end)
B(2:end) = beta + alpha*B(2:end);

% Index requested coefficients and reshape
A = reshape(A(n+1), nsize); 
B = reshape(B(n+1), nsize); 
C = reshape(C(n+1), nsize);

%% New recurrence in terms of old coefficients on standard domain
%alpha(2:end) = sqrt(b(2:end)).*c(1:end-1)./c(2:end);
%beta(2:end) = -a(1:end-1);
%gamma(3:end) = -sqrt(b(2:end-1)).*c(2:end-1)./c(1:end-2);
%alpha(1) = c(1)/sqrt(b(1));
%
%% Pluck out appropriate values of n
%alpha = alpha(n+1);
%beta = beta(n+1);
%gamma = gamma(n+1);
%
%% Assuming functions aren't scaled after mapped, this maps the recurrence
%% relation to the physical domain
%alpha = alpha./self.map_to_standard_domain.A;
%beta = (self.map_to_standard_domain.b + beta)./self.map_to_standard_domain.A;
%gamma = gamma./self.map_to_standard_domain.A;
