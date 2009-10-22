function[p] = eval_polynomial(x,alpha,beta,n,varargin);
% [P] = EVAL_POLYNOMIAL(X,ALPHA,BETA,N,{D=0,NORMALIZATION='NORMAL',SHIFT=0,SCALE=1});
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
%     SHIFT and SCALE are the affine scaling parameters. They are used for two
%     purposes:
%      - shifting X back to the standard interval
%      - normalizing the 'monic' normalization
%     
%     Monic:
%     p_{n+1} = (x-a_{n})*p_n - b_{n}*p_{n-1}
%     Normalized:
%     sqrt(b_{n+1}) p_{n+1} = (x-a_n) p_n - sqrt(b_n) p_{n-1}

global handles;
sss = handles.speclab.common.standard_scaleshift_1d.handle;
inputs = handles.common.input_schema.handle;
opt = inputs({'d','scale','shift','normalization'},{0,1,0,'normal'},[],varargin{:});

% Pre-processing:
x = x(:);
n = n(:);
d = opt.d(:);
shift = opt.shift;
scale = opt.scale;
D = max(d);
N = max(n);

x = sss(x,opt);

p = zeros([length(x) length(n) length(d)]);
p0 = zeros([length(x) N+1]);
p1 = zeros([length(x) N+1]);

if strcmpi(opt.normalization,'normal')
  beta = sqrt(beta);
  p(:,n==0,d==0) = 1/beta(1);
elseif strcmpi(opt.normalization,'monic')
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
  if strcmpi(opt.normalization,'normal')
    p1(:,qq+1) = factorial(qq)/prod(beta(1:(qq+1)));
    p1(:,qq+2) = 1/beta(qq+2)*((x-alpha(qq+1)).*p1(:,qq+1) + qq*p0(:,qq+1));
  elseif strcmpi(opt.normalization, 'monic')
    p1(:,qq+1) = factorial(qq);
    p1(:,qq+2) = (x-alpha(qq+1)).*p1(:,qq+1) + qq*p0(:,qq+1);
  end

  for q = (qq+2):N
    if strcmpi(opt.normalization,'normal')
      p1(:,q+1) = 1/beta(q+1)*((x-alpha(q)).*p1(:,q) - ...
                               beta(q)*p1(:,q-1) + ...
                               qq*p0(:,q));
    elseif strcmpi(opt.normalization,'monic')
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

p = squeeze(p);
