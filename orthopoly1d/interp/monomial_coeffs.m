function[c] = monomial_coeffs(l,W,p,k,u,f)
% monomial_coeffs -- Computes monomial coefficients for interpolation
%
% c = monomial_coeffs(l,W,p,k,u,f)
%
%     Using the 5-tuple (l,W,p,k,u) that is output from de Boor's algorithm
%     monomial_lu, this function uses the point-evaluations f to compute the
%     multimonomial coefficients necessary for interpolation. Obviously, f must
%     be an N-vector, where N is the number of rows of W (or l, or p, or u).
%
%     The interpolant may be evaluated as 
%
%       p(x) = sum_n c(n) * m(x,n),
%
%     where m(x,n) is a function that evaluates the n'th multimonomial at the
%     location x. (e.g. speclab.monomials.multimonomial is this function).
