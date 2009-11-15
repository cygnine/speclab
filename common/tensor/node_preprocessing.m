function[x,info] = node_preprocessing(x,dim)
% node_preprocessing -- array size preprocessing operators
%
% [x, info] = node_preprocessing(x, dim)
%
%     Performs some elementary preprocessing on the input x so that each column
%     corresponds to a dimension. The changes are saved in the struct output
%     'info', which can be fed into node_postprocessing to reverse the
%     operation for user output.

info.dim = dim;
info.xsize = size(x);
info.action = '';

if dim==1
  if info.xsize(2)==1  % It's already a column vector...just pass through
    % do nothing, pass
  else
    x = x(:);
    info.action = '1dresize';
  end

elseif info.xsize(2) ~= dim

  assert(info.xsize(1)==dim, 'Error: input node dimensions not compatible');

  x = x.';
  info.action = 'transpose';

else  % do nothing
end
