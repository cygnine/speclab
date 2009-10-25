function[p] = eval_polynomial_standard(x,alpha,beta,n);
% [p] = eval_polynomial_standard(x,alpha,beta,n)
%
%     Evaluates the L^2-normalized orthogonal polynomials defined by the recurrence
%     coefficients alpha and beta. Assumes alpha and beta are long enough as
%     necessary to evaluate the first max(N)'th polynomials. 
%
%     The output P is of size [length(X), length(N)].
%
%     This function is a quick, optimized version of eval_polynomial and should
%     be used instead if multiple calls to eval_polynomial are needed. The
%     functionality here is very restricted, but it's faster.
%
%     Normalized:
%     sqrt(b_{n+1}) p_{n+1} = (x-a_n) p_n - sqrt(b_n) p_{n-1}

global packages;
%sss = packages.speclab.common.standard_scaleshift_1d.handle;
%inputs = packages.labtools.input_schema.handle;
%opt = inputs({'d','scale','shift','normalization'},{0,1,0,'normal'},[],varargin{:});

% Pre-processing:
x = x(:);
n = n(:);
%if exist('d')
%  d = d(:);
%else
%  d = 0;
%end

%D = max(d);
N = max(n);

%p = zeros([length(x) length(n) length(d)]);
p = zeros([length(x) length(n)]);
p0 = zeros([length(x) N+1]);
p1 = zeros([length(x) N+1]);

beta = sqrt(beta);
%p(:,n==0,d==0) = 1/beta(1);
p(:,n==0) = 1/beta(1);

% Trivial case
if N==0; 
  return;
end

% 3D index counter for output p
%Dcount = 1;

%for qq = 0:D
  %p1(:,1:qq) = 0;
  %p1(:,qq+1) = factorial(qq)/prod(beta(1:(qq+1)));
  p1(:,1) = 1/beta(1);
  %p1(:,qq+2) = 1/beta(qq+2)*((x-alpha(qq+1)).*p1(:,qq+1) + qq*p0(:,qq+1));
  p1(:,2) = 1/beta(2)*(x-alpha(1)).*p1(:,1);

  %for q = (qq+2):N
  for q = (2):N
      p1(:,q+1) = 1/beta(q+1)*((x-alpha(q)).*p1(:,q) - ...
                               beta(q)*p1(:,q-1));
  end

  %if any(qq==d)
    %p(:,:,Dcount) = p1(:,n+1);
  %  Dcount = Dcount + 1;
  %end
  p = p1(:,n+1);

%  p0 = p1;
%end

p = squeeze(p);
