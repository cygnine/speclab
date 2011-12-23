function[w] = composition(self, u, v)
% composition -- Composes two spectral expansions
%
% w = composition(self, u, v)
%
%     If u is a length-N vector containing expansion coefficients for a
%     polynomial, and v is a length-M vector also containing expansion
%     coefficients, w contains expansion coefficients corresponding to the
%     functional composition of u with v: i.e. w is a length (N-1)*(M-1)+1
%     vector that is the composition u \circ v.

u = u(:); v = v(:);
N = length(u);
M = length(v);

% The code below works for orthonormal expansions, so let's make things orthonormal
factors = self.norm(self.range(N));
u = u.*factors;
factors = self.norm(self.range(M));
v = v.*factors;

w = zeros([(N-1)*(M-1)+1 1]);

[a,b] = self.recurrence(self.range(N+1));
b = sqrt(b);

% Temporary storage arrays
storage_n1 = w; % for stage n-1
storage_n = w;  % for stage n
storage = w;    % for current stage (n+1)

% We will do this by progressively determining the expansion of p_n \circ v,
% and then we'll iterate on n, and then add these components together at the
% end.

% for p_0 \circ v:
storage_n1(1) = 1;
% Update w:
w = w + u(1)*storage_n1;
if N==1
  return
end

% for p_1 \circ v:
storage_n(1) = -a(1);
storage_n(1:M) = storage_n(1:M) + 1/b(1)*v;
storage_n = storage_n/b(2);

% Update w:
w = w + u(2)*storage_n;
if N==2
  return
end

% Now iterate on n:
for n = 2:(N-1);
  % for p_n \circ v:
  storage(1:(1+n*(M-1))) = self.convolution(storage_n(1:((n-1)*(M-1)+1)), v);
  storage = storage - a(n)*storage_n - b(n)*storage_n1;
  storage = storage/b(n+1);

  storage_n1 = storage_n;
  storage_n = storage;
  
  % Update w:
  w = w + u(n+1)*storage;
end

factors = self.norm(self.range(length(w)));
w = w./factors;
