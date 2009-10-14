function[nodes] = canonical_nodes(N)
% [nodes] = canonical_nodes(n)
%
%     Returns the N-point canonical Fourier equidistant nodal set over the
%     interval [0,2*pi].

nodes = linspace(0,2*pi,N+1);
nodes = nodes(1:N).';