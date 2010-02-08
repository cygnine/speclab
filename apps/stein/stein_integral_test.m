function[I] = stein_integral_test(f, grid)
% stein_integral_d2_driver -- Workhorse for computing the Stein integral
%
% I = stein_integral_d2_driver(f, grid)
%
%     Assumes a gauss_grid.

%N = length(grid.theta);
%s_size = size(grid.global_s_nodes);
%contributions = zeros([N 1]);

% Evaluate all values of f:
%fx = f([grid.theta(:); grid.global_s_nodes(:)]);
%fs = reshape(fx(N+1:end), s_size);
%fx(N+1:end) = [];

cell_max = 20;
N_cells = size(grid.theta_grid,2);
N_iterations = ceil(N_cells/cell_max);

theta_contributions = zeros(size(grid.theta_grid));

N_cols= grid.N_p*cell_max;
s_grid = repmat(grid.s, [1 N_cols]);

for q = 1:N_iterations
  i1 = (q-1)*cell_max + 1;
  i2 = min(q*cell_max, N_cells);
  cell_ids = i1:i2;

  theta = grid.theta_grid(:, cell_ids);
  theta = repmat(theta(:).', [grid.N_s 1]);

  if q==N_iterations
    s_grid(:,size(theta,2)+1:end) = [];
  end

  temp = grid.w_s.'*abs(f(theta+s_grid) + f(theta-s_grid) - 2*f(theta)).^2;
  
  theta_contributions(:,i1:i2) = reshape(temp, [grid.N_p, i2-i1+1]);
end

I = grid.b*(grid.standard_element.w.'*theta_contributions)*grid.theta_scales;
