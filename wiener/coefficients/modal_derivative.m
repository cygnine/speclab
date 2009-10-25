function[inds] = modal_derivative(k,varargin)
% [inds] = modal_derivative(k,{s=1, t=0, scale=1})
%
%     Given a vector k of modal indices, returns a length(k) x 6 matrix. 
%     Each row represents the modal coefficient expansion of the derivative of
%     the Wiener rational basis function of index k. The columns of inds are
%     arranged as follows:
%
%     inds(:,1) <-----> -(k + sgn(k)) = -sgn(k)*(|k|+1)
%     inds(:,2) <-----> -k
%     inds(:,3) <-----> -(k - sgn(k)) = -sgn(k)*(|k|-1)
%     inds(:,4) <-----> k - sgn(k)    = sgn(k)*(|k|-1)
%     inds(:,5) <-----> k
%     inds(:,6) <-----> k + sgn(k)    = sgn(k)*(|k|+1)
%
%     When k=0: (N/A columns are zeroed-out)
%     inds(:,1) <-----> -1
%     inds(:,2) <-----> N/A
%     inds(:,3) <-----> N/A
%     inds(:,4) <-----> N/A
%     inds(:,5) <-----> 0
%     inds(:,6) <-----> 1
%
%     When abs(k)=1:
%     inds(:,1) <-----> -2*sgn(k)
%     inds(:,2) <-----> -sgn(k)
%     inds(:,3) <-----> N/A
%     inds(:,4) <-----> 0
%     inds(:,5) <-----> sgn(k)
%     inds(:,6) <-----> 2*sgn(k)

global packages;
wiener = packages.speclab.wiener;
opt = wiener.defaults(varargin{:});
s = opt.s;
inds = zeros([length(k) 6]);

k0 = k==0;
kn0 = not(k0);
k1 = abs(k)==1;
kk = not(or(k0,k1));

% k=0 first
inds(k0,1) = (1+sqrt(s))/sqrt(1+s);
inds(k0,5) = 2*sqrt(s-1/2);
inds(k0,6) = (1-sqrt(s))/sqrt(1+s);
inds(k0,:) = inds(k0,:)*-i/2*sqrt(s-1/2);

% all k~=0
n = abs(k(kn0))-1;
inds(kn0,1) = -s./(2*n+s+2).*(sqrt((n+2).*(n+s)) - sqrt((n+1).*(n+s+1))) + ...
                sign(k(kn0)).*(sqrt((n+1).*(n+2)) - sqrt((n+s).*(n+s+1)));
inds(kn0,6) = -s./(2*n+s+2).*(sqrt((n+2).*(n+s)) + sqrt((n+1).*(n+s+1))) + ...
                sign(k(kn0)).*(sqrt((n+1).*(n+2)) + sqrt((n+s).*(n+s+1)));
inds(kn0,1) = inds(kn0,1).*i/4.*sqrt(1-(s*(s-2))./((2*n+s+1).*(2*n+s+3)));
inds(kn0,6) = inds(kn0,6).*i/4.*sqrt(1-(s*(s-2))./((2*n+s+1).*(2*n+s+3)));
inds(kn0,5) = i*sign(k(kn0)).*sqrt((n+1).*(n+s)) - ...
             i*s*(s-1)^2./(2*(2*n+s).*(2*n+s+2)) - i*s/2;
inds(kn0,2) = i*s*(s-1)./(2*(2*n+s).*(2*n+s+2));

% General k != 0,1,-1
n = abs(k(kk))-1;
inds(kk,3) = sign(k(kk)).*(sqrt((n+s-1).*(n+s)) - sqrt((n+1).*n)) + ...
              -s./(2*n+s).*(sqrt((n+1).*(n+s-1)) - sqrt(n.*(n+s)));
inds(kk,4) = sign(k(kk)).*(sqrt((n+s-1).*(n+s)) + sqrt((n+1).*n)) + ...
              -s./(2*n+s).*(sqrt((n+1).*(n+s-1)) + sqrt(n.*(n+s)));
inds(kk,3) = inds(kk,3).*i/4.*sqrt(1-(s*(s-2))./((2*n+s+1).*(2*n+s-1)));
inds(kk,4) = inds(kk,4).*i/4.*sqrt(1-(s*(s-2))./((2*n+s+1).*(2*n+s-1)));

% special case for k=1
inds(k1,4) = i/2*sqrt((2*s-1)/(2*s+2))*(sign(k(k1))*sqrt(s)-1);
