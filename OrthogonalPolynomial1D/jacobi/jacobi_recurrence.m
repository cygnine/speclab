function[as,bs] = jacobi_recurrence(N,alpha,beta,shift,scale)

% [x,w] = jacobi_gq(N,alpha,beta)
% Returns the recurrence coefficients to generate the first N Jacobi
% polynomials defined by (alpha, beta);
%
% 20080623: acn

% Sets default values for alpha, beta
jacobi_parameters;

as = (beta^2-alpha^2)*ones([N 1]);
bs = ones([N 1]);

% Initial conditions:
as(1) = (beta-alpha)/(alpha+beta+2);
bs(1) = 2^(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+2);

for q = 2:N;
  k = q-1;
  as(q) = as(q)./((2*k+alpha+beta)*(2*k+alpha+beta+2));
  if k==1;
    bs(q) = 4*k*(k+alpha)*(k+beta)/((2*k+alpha+beta)^2*(2*k+alpha+beta+1)); 
  else
    bs(q) = 4*k*(k+alpha)*(k+beta)*(k+alpha+beta)/((2*k+alpha+beta)^2*(2*k+alpha+beta+1)*(2*k+alpha+beta-1));
  end
end

% Scale and shift:
recurrence_scaleshift;
