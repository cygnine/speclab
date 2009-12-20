function[fz] = vandermonde_interpolant(theta,f,z);
% [fz] = vandermonde_interpolant(theta,f,z);
%
%     Constructs a Fourier modal interpolant of the datapoints (theta,f) assumed
%     to be on the interval [0,2*pi] at the locations Z. 

persistent integer_range nonuniform_nodes_to_modes
if isempty(integer_range)
  from speclab.fourier import integer_range nonuniform_nodes_to_modes
end

% Force column vector
theta = theta(:);
z = z(:);
ks = integer_range(length(theta)).';

modes = nonuniform_nodes_to_modes(theta,f);

vandermonde = 1/sqrt(2*pi)*exp(i*z*ks);

fz = vandermonde*modes;
