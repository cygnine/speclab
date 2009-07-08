function[p] = eval_opoly(x,alpha,beta,n);

% function[p] = eval_opoly(x,alpha,beta,n);
% Evaluates the monic orthogonal polynomials defined by the recurrence
% coefficients alpha and beta. Assumes alpha and beta are long enough as
% necessary to evaluate the max(n)'th polynomial. 
% 
% Supports vectorization in n and x.

% 20080522: acn

% p_{n+1} = (x-a_{n})*p_n - b_{n}*p_{n-1}

% Pre-processing:
x = x(:);
n = n(:);
N = max(n);

p = zeros([length(x) N+1]);

p(:,1) = 1;
if N==0; 
  p = p(:,n+1);
  return;
end

p(:,2) = p(:,1).*(x-alpha(1));

for q=1:N;
  p(:,q+2) = (x-alpha(q+1)).*p(:,q+1) - beta(q+1)*p(:,q);
end

p = p(:,n+1);
