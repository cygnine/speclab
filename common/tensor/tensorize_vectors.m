function[x] = tensorize_vectors(varargin)
% tensorize_vectors -- transforms many 1D arrays into a grid in N-D
%
% x = tensorize_vectors(x1,x2,x3,...)
%
%     Given N vectors x1, x2, x3, ..., this function 'tensorizes' them by
%     building a grid in N-D space. The output x has N columns, and length(x1) x
%     length(x2) x length(x3) x ... rows.
%
%      Basically, this is just a wrapper for ndgrid.

N = nargin;
Ns = zeros([N 1]);

% force column vectors
for q = 1:N
  x = varargin{q};
  varargin{q} = x(:);  
  Ns(q) = length(varargin{q});
end

xs = cell([N 1]);
[xs{:}] = ndgrid(varargin{:});

x = zeros([prod(Ns) N]);
for q = 1:N
  x(:,q) = xs{q}(:);
end
