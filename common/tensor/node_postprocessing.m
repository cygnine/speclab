function[x,info] = node_postprocessing(x,dim)
% node_postprocessing -- array size preprocessing operators
%
% x = node_postprocessing(x, info)
%
%     Reverses the operation performed by node_preprocessing.The input 'info' is
%     the output from that function.

switch info.action
case ''
  % Do nothing
case '1dresize'
  x = reshape(x, info.xsize);
case 'transpose'
  x = x.';
otherwise
  error('Unrecognized postprocessing action')
end
