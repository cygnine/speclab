function[fz] = evaluate(mc,z);
% [Y] = EVALUATE(MC,Z);
%
%     For given monomial coefficients in the rows of MC, returns the evaluation
%     of that polynomial at the locations Z. Both MC and Z should be column
%     vectors. 
%
%     If MC and Z are matrices, this function is vectorized so that the same
%     operation is performed on each column of MC and Z.
%
%     Why does this even exist given Matlab's POLYVAL? Because POLYVAL isn't
%     vectorized. NB the ordering of the coefficients is reversed from POLYVAL. 

assert(size(mc,2)==size(z,2), ...
    'Error: you must give me inputs with the same number of columns');

[n,C] = size(mc);
[Z,C] = size(z);

if C==1  % Why not use polyval?

  fz = polyval(flipud(mc),z);
  return

else     % pfft, fine: the hard way

  fz = repmat(mc(1,:),[Z,1]);
  for q = 2:n
    fz = fz + z.^(q-1)*spdiags(mc(q,:),0,C,C);
  end

  return

end
