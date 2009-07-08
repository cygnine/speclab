function[nus] = calc_jacobi_moments(alpha,beta,N,gamma,delta);

% nus = calc_jacobi_moments(alpha,beta,N,gamma,delta)
% Uses Gauss-Jacobi quadrature to evaluate the first N polynomial moments of the
% Jacobi weight function (alpha,beta) against the orthogonal polynomial defined by
% (gamma,delta).
%
% 20080521: acn

[r,w] = jacobi_gq(alpha,beta,ceil((N+1)/2));
rpolys = eval_jacobi(r,gamma,delta,0:(N-1));

nus = w.'*rpolys;
