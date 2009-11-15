function[f] = newton_derivative_evaluate(x,c,varargin)
% [f] = newton_derivative_evaluate(x,c,{z=NaN, d=1})
%
%     If z is not given:
%     For nodal locations x and modal coefficients c, calculates the
%     value of the derivative of the Newton interpolant at x(1). The Horner
%     decomposition of the interpolant makes this task very easy. 
%     If c has multiple columns (say nc of them), we assume that there are multiple
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
%
%     In either case, the optional argument d determines how many derivatives to
%     take.

persistent input_schema newton_to_monomial monomial_derivative evaluate
if isempty(input_schema)
  from labtools import input_schema;
  from speclab.newton_polynomials import newton_to_monomial;
  from speclab.monomials import monomial_derivative evaluate
end
opt = input_schema({'z', 'd'}, {NaN, 1}, [], varargin{:});

[n,C] = size(c);
if and(n==1,C>1)  % I don't think you're calling this for derivatives of
  n = C;          % constants; I'm flipping some vectors around.
  C = 1;
  x = x';
  y = y';
end

%if length(varargin)<1
if isnan(opt.z);

  f = c(n,:);
  z = x(1,:);
  if opt.d>1
    error('Not coded yet')
  end
  for q = (n-1):(-1):2
    f = f.*(z-x(q,:)) + c(q,:);
  end

  return

% Need to evaluate at points inside the cells ... gets more complicated
else
  %z = varargin{1};
  z = opt.z;
  if length(z)==0
    f = [];
    return
  end

  [Z,C] = size(z);

  % We're going to cheat and translate to monomials. In the long run, this could
  % be problematic for both speed and numerical stability. But for now I'm lazy.
  monomial_coefficients = newton_to_monomial(c,x);
  
  % Find the derivative
  for q = 1:opt.d
    monomial_coefficients = monomial_derivative(monomial_coefficients);
  end

  temp = size(z);
  % Evaluate at desired points
  %f = evaluate(monomial_coefficients,z(:));
  f = evaluate(monomial_coefficients,z);
  f = reshape(f,temp);

  return
end
