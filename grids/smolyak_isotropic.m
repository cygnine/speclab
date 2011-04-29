function[x,w] = smolyak_isotropic(quadrule, q, d)
% smolyak_isotropic -- Computes a Smolyak sparse-grid quadrature rule
%
% [x,w] = smolyak_isotropic(quadrule, q, d)
%
%     Computes the d-variate `isotropic' Smolyak sparse-grid quadrature rule
%     corresponding to the family of univariate quadrature rules specified by
%     'quadrule'.
%
%     'quadrule' is a function handle with syntax quadrule(i), which returns
%     the qaudrature rule [x1,w1] at level i = 1, 2, ....
%
%     q is the level of the d-variate quadrature rule, it must be greater than
%     or equal to d. If this condition fails, the outputs are set to empty
%     arrays.
%
%     This returns the quadrature rule A(q,d), which is standard notation in the
%     sparse grid literature; x are the nodal locations and w are the weights.
%     The rule is `isotropic', meaning that the generating univariate quadrature
%     rule is identical regardless of dimension. This routine performs
%     absolutely no checking for nested quadrature rules.

persistent smolyak_indices dim subdim
if isempty(dim)
  from speclab.grids import smolyak_indices
  from speclab.common.tensor import space_dimension as dim
  from speclab.common.tensor import subspace_dimension as subdim
end

if q<d
  x = zeros([0 d]);
  w = [];
  return
end


max_index = q-d+1;

% The power of an isotropic rule: we can just pre-generate all the univariate
% rules and all that will be left is to combine them per Smolyak's construction
quad_nodes = cell([max_index 1]);
quad_weights = cell([max_index 1]);
for level=1:max_index
  [quad_nodes{level},quad_weights{level}] = quadrule(level);
end

% Now that we know the size of each quadrature rule, can allocate:
level_nodes = zeros([d 1]);
level_indices = cell([d 1]);
for level=(q-d+1):q
  %level_nodes(level-q+d) = subdim(d,level);

  level_indices{level-q+d} = smolyak_indices(d, level);
  level_indices{level-q+d} = smolyak_indices(d, level);
  indices = level_indices{level-q+d};
  for index = 1:size(indices,1)
    n_nodes = 1;
    for dimension=1:d
      n_nodes = n_nodes * length(quad_weights{indices(index,dimension)});
    end
    level_nodes(level-q+d) = level_nodes(level-q+d) + n_nodes;
  end
end
total_nodes = sum(level_nodes);
x = zeros([total_nodes d]);
w = zeros([total_nodes 1]);

temp_nodes = cell([d 1]);
temp_weights = cell([d 1]);

current_node = 1;
for level=(q-d+1):q

  %indices = smolyak_indices(d,level);
  indices = level_indices{level-q+d};

  % I see no smart way to do this other than just for-looping
  for index = 1:size(indices,1);

    % First generate this tensor-product grid
    for dimension=1:d
      temp_nodes{dimension} = quad_nodes{indices(index,dimension)};
      temp_weights{dimension} = quad_weights{indices(index,dimension)};
    end

    if d>1
      [temp_nodes{:}] = ndgrid(temp_nodes{:}); % tensor product
      [temp_weights{:}] = ndgrid(temp_weights{:}); % tensor product
    end
    rows = current_node:(current_node+numel(temp_nodes{1})-1);
    %rows = current_node:(current_node+size(temp_nodes{1},1)-1);

    w(rows) = 1;
    for dimension=1:d
      x(rows,dimension) = temp_nodes{dimension}(:);
      w(rows) = w(rows).*temp_weights{dimension}(:);
    end
    w(rows) = w(rows)*(-1)^(q-level)*nchoosek(d-1, q-level);

    current_node = rows(end)+1;
  end
end
