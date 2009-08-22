function[scale] = resolution_scaling(L,nodes,varargin)
% [scale] = resolution_scaling(L,nodes,{resolution_fraction=1})
%
%     Returns the affine scaling parameter required so that
%     ceil(resolution_fraction*lenght(nodes)) of the nodes lie inside the
%     interval [-L,L]. Mostly used for infinite-interval expansions.

global handles;
opt = handles.common.InputSchema({'resolution_fraction'}, {1}, [], varargin{:});

N = length(nodes);
delta = min([1, opt.resolution_fraction]);

Nfrac = ceil(N*delta);

if delta==1
  scale = max(abs(nodes));
  return
end

Lbig = 1e6;
Lsmall = Lbig;
% fraction inside interval [-Lnext, Lnext]
Nnodes = sum(abs(nodes*Lsmall)<=L);
while Nnodes < Nfrac
  Lbig = Lsmall;
  Lsmall = Lsmall/2;
  Nnodes = sum(abs(nodes*Lsmall)<=L);
end

% Now [Lsmall, Lbig] contains delta. Use bisection
Lmiddle = mean([Lsmall, Lbig]);
Lsep = Lbig-Lmiddle;
Ltol = 1e-6;
while Lsep>Ltol
  Nnodes = sum(abs(nodes*Lmiddle)<=L);
  if Nnodes > Nfrac
    Lsmall = Lmiddle;
  else
    Lbig = Lmiddle;
  end
  Lmiddle = mean([Lsmall, Lbig]);
  Lsep = Lbig-Lmiddle;
end

scale = Lsmall;
