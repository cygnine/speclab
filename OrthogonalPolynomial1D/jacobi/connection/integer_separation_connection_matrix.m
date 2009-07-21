function[C] = integer_separation_connection_matrix(N,alpha,beta,A,B,varargin)
% [C] = INTEGER_SEPARATION_CONNECTION_MATRIX(N,ALPHA,BETA,A,B,{NORMALIZATION='normal'})
%
%     Returns the sparse N x N connection matrix C that connects modal
%     coefficients from Jacobi polynomials of class (ALPHA,BETA) to Jacobi
%     polynomials of class (ALPHA+A, BETA+B) when A and B are integers. This
%     connection is sparse. Therefore, the matrix C is a sparse-matrix datatype. 
%
%     The optional input NORMALIZATION defines the normalization of the
%     polynomials P, which affects the matrix C.

global handles;
opt = handles.common.InputSchema({'normalization'}, {'normal'}, [],varargin{:});

coeffs = handles.speclab.OrthogonalPolynomial1D.jacobi.coefficients;

ab_min = min([A,B]);
C = spalloc(N,N,(A+B)*N - (A+B)*(A+B+1)/2);
C = C + speye(N);

ns = (0:N).';

% Explicity computing (1-r^2) x p connections is more stable than sequentially
% doing (1-r) x p and (1+r) x p.
for q = 1:ab_min
  mu = coeffs.one_minus_r_squared_times_p(ns,alpha+q,beta+q,opt);
  C_temp = spdiags(mu,[0,-1,-2],N,N).';
  C = C_temp*C;
end

alpha = alpha + ab_min;
beta = beta + ab_min;

% Compute remaining (1 \pm r) x p connections
if A>ab_min
  Q = A - ab_min;
  for q = 1:Q
    mu = coeffs.one_minus_r_times_p(ns,alpha+q,beta,opt);
    C_temp = spdiags(mu,[0,-1],N,N).';
    C = C_temp*C;
  end
elseif B>ab_min
  Q = B - ab_min;
  for q = 1:Q
    mu = coeffs.one_plus_r_times_p(ns,alpha,beta+q,opt);
    C_temp = spdiags(mu,[0,-1],N,N).';
    C = C_temp*C;
  end
end
