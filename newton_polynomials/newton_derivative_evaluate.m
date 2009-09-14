function[f] = newton_derivative_evaluate(x,c,varargin)
% [F] = NEWTON_DERIVATIVE_EVALUATE(X,C)
% [F] = NEWTON_DERIVATIVE_EVALUATE(X,C,Z)
%
%     If Z is not given:
%     For nodal locations X and modal coefficients C, calculates the
%     value of the derivative of the Newton interpolant at X(1). The Horner
%     decomposition of the interpolant makes this task very easy. 
%     If C has multiple columns (say NC of them), we assume that there are multiple
%     interpolants whose derivatives must be evaluated simultaneously.
%
%     For vectors:
%     length(C) = n
%     length(X) = n-1  (any additional nodal locations are ignored)
%
%     If Z is given:
%     Computes the derivatives of the functions at the locations Z. 
%     [Z,C] = size(z);
%     Z is just the number of evaluations per cell. Unless this is going to be
%     constant across cells, you might be better off feeding one cell at a time
%     to the this function and using a parfor to parallelize the outer loop. 

global handles;
newton = handles.speclab.newton_polynomials;
mono = handles.speclab.monomials;

[n,C] = size(c);
if and(n==1,C>1)  % I don't think you're calling this for derivatives of
  n = C;          % constants; I'm flipping some vectors around.
  C = 1;
  x = x';
  y = y';
end

if length(varargin)<1

  f = c(n,:);
  z = x(1,:);
  for q = (n-1):(-1):2
    f = f.*(z-x(q,:)) + c(q,:);
  end

  return

% Need to evaluate at points inside the cells ... gets more complicated
else
  z = varargin{1};
  if length(z)==0
    f = [];
    return
  end

  [Z,C] = size(z);

  % We're going to cheat and translate to monomials. In the long run, this could
  % be problematic for both speed and numerical stability. But for now I'm lazy.
  monomial_coefficients = newton.newton_to_monomial(c,x);
  
  % Find the derivative
  monomial_coefficients = mono.monomial_derivative(monomial_coefficients);

  temp = size(z);
  % Evaluate at desired points
  f = mono.evaluate(monomial_coefficients,z(:));
  f = reshape(f,temp);

  return
end
