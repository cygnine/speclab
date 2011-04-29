function[z, x] = rkhs_lattice1_generating_vector(n, d, varargin)
% rhks_lattice1_generating_vector -- generating vector for rank-1 lattice rule
%
% z = rkhs_lattice1_generating_vector(n, d, {...})
%
%     Calculates the generating vector for the rank-1 lattice rules for
%     worst-case error in tensor-product reproducing kernel Hilbert spaces. This
%     function assumes that n is prime and no checking is performed to ensure
%     this.
%
%     As of now, this is a direct approach, consuming O(d*n^2) time. 
%     TODO: implement the O(d*n*log(n)) time approach.

persistent strict_inputs
if isempty(strict_inputs)
  from labtools import strict_inputs
end

gamma_default = 0.9.^(1:d);
opt = strict_inputs({'alpha', 'gamma', 'compute_nodes'}, {2, gamma_default, false}, [], varargin{:});

omega = generate_omega(opt.alpha);

% storage for the omega factors
temp_storage = zeros([n 1]);
% storage for discrepancy values
e_storage = ones([n 1]);
e_values = ones([1 n-1]);

k = (0:n-1).';
x_temp = k/n;
if opt.compute_nodes
  x = zeros([n d]);
  x(:,1) = x_temp;
else
  x = [];
end

prefer_speed = n<1000;

% Initial condition
z = zeros([d 1]);
z(1) = 1;

if prefer_speed % We can sacrifice storage for speed
  ks_matrix = omega(mod(k*k(2:end).'/n, 1));
else % We probably need to sacrifice speed for memory issues
end

for dim = 2:d
  e_storage = e_storage.*(1+opt.gamma(dim-1)*omega(x_temp));

  if prefer_speed % The storage version -- easy to implement
    e_values = e_storage'*(1 + opt.gamma(dim)*ks_matrix)/n - 1;
    [garbage, min_index] = min(e_values);
    x_temp = mod(k*k(min_index+1)/n,1);
    z(dim) = min_index;
  else % The non-storage version -- a little harder
    prev_evalue = Inf;
    for qq = 1:(n-1);
      y_temp = mod(k*qq/n,1);
      e_values(qq) = -1 + 1/n*e_storage'*(1+opt.gamma(dim)*omega(y_temp));
      if e_values(qq) < prev_evalue
        x_temp = y_temp;
        z(dim) = qq;
        prev_evalue = e_values(qq);
      end
    end
  end
  
  if opt.compute_nodes
    x(:,dim) = x_temp;
  end
end

%%%%%%%%%%%%%%% SUBFUNCTION %%%%%%%%%%%%%%% 
function[omega] = generate_omega(alpha)
% generate_omega -- Function handle for evaluating the `variable' part of the kernel in a Korobov space
%
% omega = omega(alpha)

persistent bernoulli
if isempty(bernoulli)
  from labtools.specfun import bernoulli_polynomial as bernoulli
end

if mod(alpha,2)<1e-10
  alpha = round(alpha);
  omega_sign = (-1)^(1+mod(alpha/2,2));
  afact = factorial(alpha);
  omega = @(x) (2*pi)^alpha*omega_sign/afact*bernoulli(x, alpha);
else
  error('Not yet implemented')
end
