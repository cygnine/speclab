function[output_modes] = modal_padding(modes,N)
% [OUTPUT_MODES] = MODAL_PADDING(MODES,N)
%
%     Symmetrically zero-pads the modal vector MODES so that it has length N.
%     This is basically a utility function that can be used to take a low-order
%     approximation and use the FFT to interpolate it at a high-order equispaced
%     collection of nodes.

N_coarse = length(modes);
% Force column vector
modes = modes(:);

if N_coarse >= N
  fprintf(['Cannot zero-pad...you asked for a shorter vector', ...
  'than already present\n']);
  NLeft = 0;
  NRight = 0;

elseif mod(N-N_coarse,2)==0

  NLeft = (N-N_coarse)/2;
  NRight = NLeft;

elseif mod(N_coarse,2)==0
  
  NLeft = (N-N_coarse-1)/2;
  NRight = NLeft+1;

else

  NLeft = (N-N_coarse+1)/2;
  NRight = NLeft-1

end

output_modes = [zeros([NLeft,1]); modes; zeros([NRight,1])];
