function[S] = wiener_stiffness_matrix(N,varargin)
% [S] = wiener_stiffness_matrix(N,{s=1, t=0, scale=1})
%
%     Returns the N x N stiffness matrix for the Wiener rational functions. If N
%     is odd, the indexing is negatively-biased. S is sparse, and contains ~3*N
%     entries if s=1, and ~6*N entries otherwise.
%
%     TODO: this is hardly the most efficient way to construct this....

global handles;
wiener = handles.speclab.wiener;
opt = wiener.defaults(varargin{:});

ks = handles.speclab.common.integer_range(N+4);

if opt.s==1
  S = spalloc(N+4,N+4,3*(N+4));
else
  S = spalloc(N+4,N+4,6*(N+4));
end
inds = wiener.coefficients.modal_derivative(ks,opt);
inds = inds/opt.scale;

zeroind = find(ks==0);

for q = 3:(N+2)
  k = ks(q);
  if k==0
    S([q-1,q,q+1],q) = inds(q,[1,5,6]).';
  elseif abs(k)==1
    cinds = [zeroind-2*sign(k), zeroind-sign(k), zeroind, ...
             zeroind+sign(k), zeroind+2*sign(k)];
    S(cinds,q) = inds(q,[1,2,4,5,6]).';
  else
    cinds = zeroind + [-(k+sign(k)), -k, -(k-sign(k)), ...
                       k-sign(k), k, k+sign(k)];
    S(cinds,q) = inds(q,:).';
  end
end

S = S(3:(N+2),3:(N+2));
