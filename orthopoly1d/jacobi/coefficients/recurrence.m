function[a,b] = recurrence(N,varargin)
% [a,b] = recurrence(N,{alpha=-1/2,beta=-1/2})
%
%     Returns the recurrence coefficients to generate the first N Jacobi
%     polynomials defined by (alpha, beta). 

% Sets default values for alpha, beta
global handles;
opt = handles.speclab.orthopoly1d.jacobi.defaults(varargin{:});
[alpha,beta] = deal(opt.alpha,opt.beta);

N = double(N);
a = (beta^2-alpha^2)*ones([N 1]);
b = ones([N 1]);

% Initial conditions:
a(1) = (beta-alpha)/(alpha+beta+2);
b(1) = 2^(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+2);

for q = 2:N;
  k = q-1;
  a(q) = a(q)./((2*k+alpha+beta)*(2*k+alpha+beta+2));
  if k==1;
    b(q) = 4*k*(k+alpha)*(k+beta)/((2*k+alpha+beta)^2*(2*k+alpha+beta+1)); 
  else
    b(q) = 4*k*(k+alpha)*(k+beta)*(k+alpha+beta)/((2*k+alpha+beta)^2*(2*k+alpha+beta+1)*(2*k+alpha+beta-1));
  end
end
