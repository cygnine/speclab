function[scale,shift] = box_to_scaleshift(x_box)
% box_scaleshift -- Returns affine scale, shift parameters for a multidimensional box
%
% [scale,shift] = box_to_scaleshift(x_box)
%
%     Returns the affine scale and shift parameters for the dim-dimensional box
%     specified by the (2 x dim) or (dim x 2) array x_box.
%
%     If dim==2, the dimensions of x_box are assumed (2 x dim).
%
%     Example:
%
%        [scale, shift] = box_to_scaleshift([-1   0; 
%                                             1   pi]);
% 
%     returns scale = [1 pi/2], shift = [0 pi/2]

box_size = size(x_box);

if box_size(1)==1
  assert(box_size(2)==2, 'One of the box dimensions must be 2');

  shift = mean(x_box);
  scale = diff(x_box)/2;
elseif box_size(2)==1
  assert(box_size(1)==2, 'One of the box dimensions must be 2');

  shift = mean(x_box);
  scale = diff(x_box)/2;
else

  if box_size(2)==2

    scale = diff(x_box, 1, 1)/2;
    shift = mean(x_box, 1);

  elseif box_size(1)==2

    scale = diff(x_box, 1, 2).'/2;
    shift = mean(x_box,2).';

  else

    error('One of the input dimensions must be 2');
  end
end
