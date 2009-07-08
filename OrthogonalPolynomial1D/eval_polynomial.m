function[p] = eval_polynomial(x,alpha,beta,n,varargin);
% [P] = EVAL_POLYNOMIAL(X,ALPHA,BETA,N,{D:0,NORMALIZATION='NORMAL'});
%
%     Evaluates the normalized orthogonal polynomials defined by the recurrence
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

global handles;
sss = handles.speclab.common.standard_scaleshift_1d;
opt = handles.common.InputSchema({'d','scale','shift','normalization'},{0,1,0,'normal'},[],varargin{:});

% Pre-processing:
x = x(:);
n = n(:);
d = opt.d(:);
shift = opt.shift;
scale = opt.scale;
D = max(d);
N = max(n);

x = sss(x,scale,shift);

p = zeros([length(x) length(n) length(d)]);
p0 = zeros([length(x) N+1]);
p1 = zeros([length(x) N+1]);

switch opt.normalization
case 'normal'
  beta = sqrt(beta);
  p(:,n==0,d==0) = 1/beta(1);
case 'monic'
  p(:,n==0,d==0) = 1;
end

% Trivial case
if N==0; 
  return;
end

% 3D index counter for output p
Dcount = 1;

for qq = 0:D
  p1(:,1:qq) = 0;
  switch opt.normalization
  case 'normal'
    p1(:,qq+1) = factorial(qq)/prod(beta(1:(qq+1)));
    p1(:,qq+2) = 1/beta(qq+2)*((x-alpha(qq+1)).*p1(:,qq+1) + qq*p0(:,qq+1));
  case 'monic'
    p1(:,qq+1) = factorial(qq);
    p1(:,qq+2) = (x-alpha(qq+1)).*p1(:,qq+1) + qq*p0(:,qq+1);
  end

  for q = (qq+2):N
    switch opt.normalization
    case 'normal'
      p1(:,q+1) = 1/beta(q+1)*((x-alpha(q)).*p1(:,q) - ...
                               beta(q)*p1(:,q-1) + ...
                               qq*p0(:,q));
    case 'monic'
      p1(:,q+1) = ((x-alpha(q)).*p1(:,q) - ...
                               beta(q)*p1(:,q-1) + ...
                               qq*p0(:,q));
    end
  end

  if any(qq==d)
    p(:,:,Dcount) = p1(:,n+1)/(scale^qq);
    Dcount = Dcount + 1;
  end

  p0 = p1;
end

p = squeeze(p)/sqrt(scale);
