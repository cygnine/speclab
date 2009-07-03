function[fz] = vandermonde_interpolant(theta,f,z);
% [FZ] = VANDERMONDE_INTERPOLANT(THETA,F,Z);
%
%     Constructs a Fourier modal interpolant of the datapoints (theta,f) assumed
%     to be on the interval [0,2*pi]. 

global handles;
fourier = handles.bases.fourier;

% Force column vector
theta = theta(:);
z = z(:);
ks = fourier.integer_range(length(theta)).';

modes = fourier.nonuniform_nodes_to_modes(theta,f);

vandermonde = 1/sqrt(2*pi)*exp(i*z*ks);

fz = vandermonde*modes;
