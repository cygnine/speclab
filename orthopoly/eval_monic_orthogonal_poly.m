function[p] = eval_monic_orthogonal_poly(x,alpha,beta,n);
% [P] = EVAL_NORMALIZED_ORTHOGONAL_POLY(X,ALPHA,BETA,N,{D:0,NORMALIZATION='NORMAL'});
%
%     Evaluates the monic orthogonal polynomials defined by the recurrence
%     coefficients ALPHA and BETA. Assumes ALPHA and BETA are long enough as
%     necessary to evaluate the first max(N)'th polynomials. 
%
%     The output P is of size [length(X), length(N)].
%
%     The optional input D is a vector of whole numbers designating which
%     derivatives are to be evaluated. The default 0 indicates just to evaluate
%     the polynomials. Vectorization in D is realized via the third dimension in
%     the output P.
%
%     The second optional input NORMALIZATION determines which normalization to
%     use. The possible values are:
%
%       'normal':  L^2(w) normalized polynomials, where w is the orthogonal
%         weight function
%     
%     Monic:
%     p_{n+1} = (x-a_{n})*p_n - b_{n}*p_{n-1}
%     Normalized:
%     sqrt(b_{n+1}) p_{n+1} = (x-a_n) p_n - sqrt(b_n) p_{n-1}

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
