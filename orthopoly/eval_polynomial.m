function[p] = eval_polynomial(x,alpha,beta,n,varargin);
% [p] = eval_polynomial(x,alpha,beta,n,{d=0,normalization='normal',shift=0,scale=1});
%
%     Evaluates the normalized orthogonal polynomials defined by the recurrence
%     coefficients ALPHA and BETA. Assumes ALPHA and BETA are long enough as
%     necessary to evaluate the first max(N)'th polynomials. 
%
%     The output p is of size [length(x), length(n)].
%
%     The optional input d is a vector of whole numbers designating which
%     derivatives are to be evaluated. The default 0 indicates just to evaluate
%     the polynomials. Vectorization in d is realized via the third dimension in
%     the output p.
%
%     The second optional input 'normalization' determines which normalization to
%     use. The possible values are:
%
%       'normal':  L^2(w) normalized polynomials, where w is the weight function
%
%     'shift' and 'scale' are the affine scaling parameters. They are used for two
%     purposes:
%      - shifting x back to the standard interval
%      - normalizing the 'monic' normalization
%     
%     Monic:
%     p_{n+1} = (x-a_{n})*p_n - b_{n}*p_{n-1}
%
%     Normalized:
%     sqrt(b_{n+1}) p_{n+1} = (x-a_n) p_n - sqrt(b_n) p_{n-1}

persistent input_schema sss
if isempty(input_schema)
  from labtools import input_schema
  from speclab.common import standard_scaleshift_1d as sss
end
opt = input_schema({'d','scale','shift','normalization'},{0,1,0,'normal'},[],varargin{:});

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

switch lower(opt.normalization);
case 'normal'
  for qq = 0:D
    p1(:,1:qq) = 0;
    p1(:,qq+1) = factorial(qq)/prod(beta(1:(qq+1)));
    p1(:,qq+2) = 1/beta(qq+2)*((x-alpha(qq+1)).*p1(:,qq+1) + qq*p0(:,qq+1));

    for q = (qq+2):N
      p1(:,q+1) = 1/beta(q+1)*((x-alpha(q)).*p1(:,q) - ...
                               beta(q)*p1(:,q-1) + ...
                               qq*p0(:,q));
    end

    if any(qq==d)
      p(:,:,Dcount) = p1(:,n+1)/(scale^qq);
      Dcount = Dcount + 1;
    end

    p0 = p1;
  end
  
case 'monic'

  for qq = 0:D
    p1(:,1:qq) = 0;
    p1(:,qq+1) = factorial(qq);
    p1(:,qq+2) = (x-alpha(qq+1)).*p1(:,qq+1) + qq*p0(:,qq+1);

    for q = (qq+2):N
      p1(:,q+1) = ((x-alpha(q)).*p1(:,q) - ...
                               beta(q)*p1(:,q-1) + ...
                               qq*p0(:,q));
    end

    if any(qq==d)
      p(:,:,Dcount) = p1(:,n+1)/(scale^qq);
      Dcount = Dcount + 1;
    end

    p0 = p1;
  end

otherwise
  error(strcat('Normalization type ''', opt.normalization, ''' not recognized'));
end

p = squeeze(p);

% for qq = 0:D
%   p1(:,1:qq) = 0;
%   if strcmpi(opt.normalization,'normal')
%     p1(:,qq+1) = factorial(qq)/prod(beta(1:(qq+1)));
%     p1(:,qq+2) = 1/beta(qq+2)*((x-alpha(qq+1)).*p1(:,qq+1) + qq*p0(:,qq+1));
%   elseif strcmpi(opt.normalization, 'monic')
%     p1(:,qq+1) = factorial(qq);
%     p1(:,qq+2) = (x-alpha(qq+1)).*p1(:,qq+1) + qq*p0(:,qq+1);
%   end
% 
%   for q = (qq+2):N
%     if strcmpi(opt.normalization,'normal')
%       p1(:,q+1) = 1/beta(q+1)*((x-alpha(q)).*p1(:,q) - ...
%                                beta(q)*p1(:,q-1) + ...
%                                qq*p0(:,q));
%     elseif strcmpi(opt.normalization,'monic')
%       p1(:,q+1) = ((x-alpha(q)).*p1(:,q) - ...
%                                beta(q)*p1(:,q-1) + ...
%                                qq*p0(:,q));
%     end
%   end
% 
%   if any(qq==d)
%     p(:,:,Dcount) = p1(:,n+1)/(scale^qq);
%     Dcount = Dcount + 1;
%   end
% 
%   p0 = p1;
% end
% 
% p = squeeze(p);
