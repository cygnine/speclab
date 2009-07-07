function[modes] = nonuniform_nodes_to_modes(theta,f);
% [MODES] = NONUNIFORM_NODES_TO_MODES(THETA,F);
%
%     Assuming that theta lives on [-pi,pi], this function uses Vandermonde
%     interpolation to compute the modal coefficients corresponding to the
%     Fourier interpolant passing through the points (THETA,F).  The modes are
%     for L^2 orthonormalized canonical Fourier basis functions. The optional
%     output KS are the integer indices corresponding to the basis functions. 

global handles;
fourier = handles.speclab.fourier;

% Force column vector
theta = theta(:);
f = f(:);

ks = fourier.integer_range(length(theta)).';

vandermonde = 1/sqrt(2*pi)*exp(i*theta*ks);

modes = inv(vandermonde)*f;
