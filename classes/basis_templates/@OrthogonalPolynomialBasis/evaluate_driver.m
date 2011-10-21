function[p] = evaluate_driver(x,alpha,beta,n,d);
% evaluate_driver -- Evaluates orthogonal polynomials
%
% p = evaluate_driver(x,alpha,beta,n,d);
%
%     Evaluates the normalized orthogonal polynomials defined by the recurrence
%     coefficients alpha and beta. Assumes alpha and beta are long enough as
%     necessary to evaluate the first max(n)'th polynomials. 
%
%     The output p is of size [length(x), length(n)].
%
%     The optional input d is a vector of whole numbers designating which
%     derivatives are to be evaluated. The default 0 indicates just to evaluate
%     the polynomials. Vectorization in d is realized via the third dimension in
%     the output p.
%
%     This function returns evaluations of the L^2-orthonormal polynomials, and
%     assumes that the inputs x are located on the standard domain of the
%     orthogonal polynomial family.
%
%     This is a static method and can be invoked without instantiating its
%     class.

% Pre-processing:
x = x(:);
n = n(:);
%d = opt.d(:);
D = max(d);
N = max(n);

p = zeros([length(x) length(n) length(d)]);
p0 = zeros([length(x) N+1]);
p1 = zeros([length(x) N+1]);

beta = sqrt(beta);
p(:,n==0,d==0) = 1/beta(1);

% Trivial case
if N==0; 
  return;
end

% 3D index counter for output p
Dcount = 1;

for qq = 0:D
  p1(:,1:qq) = 0;
  p1(:,qq+1) = factorial(qq)/prod(beta(1:(qq+1)));
  p1(:,qq+2) = 1/beta(qq+2)*((x-alpha(qq+1)).*p1(:,qq+1) + qq*p0(:,qq+1));

  for q = (qq+2):N
    p1(:,q+1) = 1/beta(q+1)*((x-alpha(q)).*p1(:,q) - ...
                             beta(q)*p1(:,q-1) + ...
                             qq*p0(:,q));
  end

  flags = find(d==qq);
  for qqq =1:length(flags)
    p(:,:,flags(qqq)) = p1(:,n+1);
    Dcount = Dcount + 1;
  end
  %if any(qq==d)
  %  p(:,:,Dcount) = p1(:,n+1);
  %  Dcount = Dcount + 1;
  %end

  p0 = p1;
end

if size(p,3)==1
  p = squeeze(p);
end
