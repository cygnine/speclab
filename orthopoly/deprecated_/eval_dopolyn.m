function[dp] = eval_dopolyn(x,alpha,beta,n);

% function[p] = eval_dopolyn(x,alpha,beta,n);
% Evaluates the derivative of the normalized orthogonal polynomials defined by the recurrence
% coefficients alpha and beta. Assumes alpha and beta are long enough as
% necessary to evaluate the max(n)'th polynomial. 
% 
% Supports vectorization in n and x.
%
% p_{n+1} = (x-a_{n})*p_n - b_{n}*p_{n-1}
% p'{n+1} = p_n + (x-a{n})*p'_n - b_{n}*p'_{n-1}

% sqrt(b_{n+1}) p_{n+1} = (x-a_n)*p_n - sqrt(b_n)*p_{n-1}
% sqrt(b_{n+1}) p'_{n+1} = p_n + (x-a_n)*p'_n - sqrt(b_n)*p'_{n-1}

% 20080523: acn

% Pre-processing:
x = x(:);
n = n(:);
N = max(n);

p = zeros([length(x) N+1]);
dp = p;

p(:,1) = 1/sqrt(beta(1));
dp(:,1) = 0;
if N==0; 
  return;
end

p(:,2) = p(:,1).*(x-alpha(1));
dp(:,2) = p(:,1);

for q=1:N;
  % Normalization of previous polynomial:
  p(:,q+1) = p(:,q+1)/sqrt(beta(q+1));
  dp(:,q+1) = dp(:,q+1)/sqrt(beta(q+1));
  % Computation of next orthogonal polynomial:
  p(:,q+2) = (x-alpha(q+1)).*p(:,q+1) - sqrt(beta(q+1))*p(:,q);
  dp(:,q+2) = p(:,q+1) + (x-alpha(q+1)).*dp(:,q+1) - sqrt(beta(q+1))*dp(:,q);
end

dp(:,end) = dp(:,end)/sqrt(beta(N+2));
dp = dp(:,n+1);
