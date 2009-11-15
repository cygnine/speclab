function[a] = linear_to_array_indexing(n,varargin)
% linear_to_array_indexing -- transforms linear indexing to array indexing
% 
% [a] = linear_to_array_indexing(n, {dim = 1})
%
%     Computes (0-based) array indexing from (0-based) linear indexing. The
%     output a is a size(n,1) x dim matrix, where each row corresponds to the
%     0-based array indexing for the given linear index n.
%
%     Examples:
%          n = [0,1,2,3,4,5,6], dim = 2 ----> [0 , 0
%                                              1 , 0
%                                              0 , 1
%                                              2 , 0
%                                              1 , 1
%                                              0 , 2
%                                              3 , 0]
%

input_schema = from_as('labtools', 'input_schema');
opt = input_schema({'dim'}, {1}, [], varargin{:});
dim = opt.dim;

n = n(:);
n = n+1;
n_max = max(n);

% Find the largest polynomial order we'll need:
N = 0;
while nchoosek(N+dim,dim)<n_max;
  N = N + 1;
end

% For right now, I don't see an easier way to do this than to just generate all
% the tuples up to max(n) and then pluck out the ones we want
a = zeros([nchoosek(N+dim,dim), dim]);

row_id = 1;
a(row_id,:) = zeros([1 dim]);
row_id = 2;

if N == 0
  a = zeros([size(n,1) dim]);
else
  for q = 1:N
    current_row = zeros([1 dim]);
    current_row(1) = q;
    a(row_id,:) = current_row;
    row_id = row_id + 1;

    % "The traveling ones-man method"
    onesman_home = 1;
    onesman_location = 1;

    finished = false;

    while not(finished);
      onesman_pilgrimage();
    end
  end

  a = a(n,:);
end

function[] = onesman_pilgrimage();
  while onesman_location < dim
    onesman_location = onesman_location + 1;
    current_row(onesman_location-1) = current_row(onesman_location-1) - 1;
    current_row(onesman_location) = current_row(onesman_location) + 1;
    a(row_id,:) = current_row;
    row_id = row_id + 1;
  end

  if onesman_home + 1 == dim 
    % Then make all the other onesman in column dim-1 travel as well
    while current_row(onesman_home)>0
      current_row(end) = current_row(end) + 1;
      current_row(end-1) = current_row(end-1) - 1;
      a(row_id,:) = current_row;
      row_id = row_id + 1;
    end
  end

  if current_row(end)==q
    finished = true;
    return % done!
  end

  % Now update new home for (next) onesman
  % There must exist a zero in some column; find the last consecutive one from
  % the right

  columns = find(current_row, 2, 'last');
  current_row(columns(1)) = current_row(columns(1)) - 1;
  current_row(columns(1)+1) = current_row(end) + 1;
  current_row(end) = 0;
  a(row_id,:) = current_row;
  row_id = row_id + 1;

  onesman_home = columns(1)+1;
  onesman_location = columns(1)+1;
end

end
