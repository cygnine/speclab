function[v] = uinv_util(u,v,d,k)
% uinv_util -- A utility to perform inversion of a matrix
%
% v = uinv_util(u,v,d,k)
%
%     This function performs v <---- inv(u)*v, except that v is a vector that
%     stores the block-diagonal elements of a matrix W. If we had the full
%     matrix W, this is simply W = inv(u)*W. However, we only care about the
%     block-diagonal elements stored in v. The scalar d is the spatial dimension
%     of the problem and the vector k indicates how blocks of the matrix W are
%     stored in the vector v.
%
%     This is a utility for least_coeffs. See least_lu for further explanation
%     of the inputs.
%
%     This function is no longer needed.

persistent subdim invu
if isempty(subdim)
  from speclab.common.tensor import subspace_dimension as subdim
  from labtools.linalg import triu_back_substitute as invu
end

Nrows = size(u,2);  % Number of rows of the W matrix
kmax = max(k);
flags = diff([0; find(diff(k)); Nrows]); % how many rows for each degree
current_row = Nrows;
v_position = length(v);

for kvalue = kmax:-1:0
  colspan = subdim(d,kvalue);
  degree_indices = (v_position-flags(kvalue+1)*colspan+1):v_position;
  degree_indices = reshape(degree_indices.', [colspan flags(kvalue+1)]).';

  row_indices = (current_row - flags(kvalue+1) + 1):current_row;

  % Now this is just a regular inversion of a submatrix of u
  v(degree_indices) = invu(u(row_indices,row_indices), v(degree_indices));
  % The above line could also just be:
  % v(degree_indices) = inv(u(row_indices,row_indices))*v(degree_indices);

  v_position = degree_indices(1,1)-1;
  current_row = row_indices(1)-1;
end
