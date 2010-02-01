function[I] = stein_integral_d2_driver(f, grid);
% stein_integral_d2_driver -- Workhorse for computing the Stein integral
%
% I = stein_integral_d2_driver(f, grid)

N = length(grid.theta);
s_size = size(grid.global_s_nodes);
contributions = zeros([N 1]);

% Evaluate all values of f:
fx = f([grid.theta(:); grid.global_s_nodes(:)]);
fs = reshape(fx(N+1:end), s_size);
fx(N+1:end) = [];

% Combine them in funky ways
for q = 1:N
  f_plus = circshift(fs, [0, 1-q]);
  f_minus = fliplr(flipud(f_plus));

  func = abs(f_plus + f_minus - 2*fx(q)).^2;

  contributions(q) = sum(sum(grid.weights.*func)); 
end

I = grid.b*sum(contributions.*grid.theta_w);
