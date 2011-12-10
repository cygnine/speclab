function[ek] = triple_product(self, N, K)
% triple_product -- weighted integral of triple product of polynomials
%
% ek = triple_product(self, N, K)
%
%     Computes the triple product coefficients
%
%       e^{(k)}_{n,m} = \langle \pi_n \pi_m, \pi_k \rangle_\omega,
%
%     for all 1 <= n,m <= N, and 0 <= k <= K.
%
%     The output ek is cell array of length K+1, where ek{p} is a sparse 2D
%     array containing the coefficients for k = p-1. (It is symmetric, but this
%     potential savings is not exploited.)

persistent spdiag
if isempty(spdiag)
  from labtools import spdiag
end

[a,b] = self.recurrence(0:N+K-1);
b = sqrt(b);

ek = cell([K+1 1]);

% k = 0 is easy:
ek{1} = 1/b(1)*speye(N+K);
if K == 0
  return
end

% k = 1 is not too hard:
p = 2;
ek{p} = spalloc(N+K,N+K,(N+K)*2);
ek{p}(1:end-1,:) = spdiag(b(2:end))*ek{p-1}(2:end,:);
ek{p}(1:end-1,:) = ek{p}(1:end-1,:) + (spdiag(a(1:end-1)) - a(p-1))*ek{p-1}(1:end-1,:);
ek{p}(2:end,:) = ek{p}(2:end,:) + spdiag(b(2:end))*ek{p-1}(1:end-1,:);
ek{p} = ek{p}/b(p);

if K == 1
  ek{1} = ek{1}(1:N,1:N);
  ek{2} = ek{2}(1:N,1:N);
  return
end

for k = 2:K
  p = k + 1;

  ek{p} = spalloc(N+K,N+K,(N+K)*p);
  ek{p}(1:end-1,:) = spdiag(b(2:end))*ek{p-1}(2:end,:);
  ek{p}(1:end-1,:) = ek{p}(1:end-1,:) + (spdiag(a(1:end-1)) - a(p-1))*ek{p-1}(1:end-1,:);
  ek{p}(2:end,:) = ek{p}(2:end,:) + spdiag(b(2:end))*ek{p-1}(1:end-1,:);
  ek{p}(1:end-1,:) = ek{p}(1:end-1,:) - b(p-1)*ek{p-2}(1:end-1,:);
  ek{p} = ek{p}/b(p);
end

for p = 1:(K+1);
  ek{p} = ek{p}(1:N,1:N);
end

% Post-processing for scaling
NK = max(N,K+1);
scales = self.scale_functions(ones([1 NK]), self.range(NK));
rowcolscale = spdiag(scales(1:N));

% Apply to each row, column, level:
for k = 0:K
  p = k+1;
  ek{p} = scales(p)*rowcolscale*ek{p}*rowcolscale;
end
