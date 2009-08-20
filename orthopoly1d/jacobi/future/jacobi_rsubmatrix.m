function[R] = jacobi_rmatrix(N,a,b,indicator,padding)

% [R] = jacobi_rmatrix(N,a,b,indicator,padding)
% Creates the R matrix which translates between the modes corresponding to (a,b)
% and those for (a-1,b-1). R is (N+padding-1) x (N+padding-1). padding is an
% optional argument, default value of 1
% indicator is a mandatory argument with the following values:
%      'a' : The R matrix demotes a
%            (1-r)*P^(a,b) = R*P^(a-1,b)
%      'b' : The R matrix demotes b
%            (1+r)*P^(a,b) = R*P^(a,b-1)
% In both cases the polynomials are orthonormal polynomials.
%
% 20080812 -- acn

if ~exist('padding', 'var');
  padding=1;
elseif padding<1;
  warning('Mininum value of padding is 1...setting padding to 1');
else
  padding = round(padding);
end

R = spalloc(N+padding-1,N+padding-1,2*(N+padding-1));

gs = jacobiconst_gammatilde(0:(N-1),a,b);

switch indicator;
case 'a';
  gs = jacobiconst_gammatilde(0:(N-1),a,b);
  gs = [gs(:,1) [0; -gs(1:(N-1),2)]];
%  for q = 1:(N-1);
%    R(q,q) = gs(q,1);
%    R(q,q+1) = -gs(q,2);
%  end
%  R(N,N) = gs(N,1);
case 'b';
  gs = jacobiconst_gammatilde(0:(N-1),b,a);
  gs = [gs(:,1) [0; gs(1:(N-1),2)]];
%  for q = 1:(N-1);
%    R(q,q) = gs(q,1);
%    R(q,q+1) = gs(q,2);
%  end
%  R(N,N) = gs(N,1);
otherwise;
  error('Bad value for indicator argument');
end

R = spdiags(gs,[0 1],N,N);
