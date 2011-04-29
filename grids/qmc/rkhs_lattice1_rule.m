function[x] = rkhs_lattice1_rule(n,d, varargin)
% rhks_lattice1_rule -- rank-1 lattice rule for a reproducing kernel Hilbert space
%
% x = rhks_lattice1_rule(n,d,{...})
%
%     Produces the rank-1 lattice rule that minimizes the discrepancy error for
%     integrands in a reproducing kernel Hilbert space. The space is defined by
%     the weights gamma, which should satisfy some convergence properties in
%     order for the rule to produce no larger than the worst-case error.
%
%     The output x is an n x d matrix: n points in d dimensions on the unit
%     hypercube [0,1]^d.

persistent generating_vector strict_inputs spdiag
if isempty(strict_inputs)
  from labtools import strict_inputs spdiag
  from speclab.grids.qmc import rkhs_lattice1_generating_vector as generating_vector
end

varargin{end+1} = 'compute_nodes';
varargin{end+1} = true;
[z, x] = generating_vector(n, d, varargin{:});

x = spdiag((0:(n-1))/n)*repmat(z(:).', [n 1]);
x = mod(x,1);
