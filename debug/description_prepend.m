function[descriptions] = description_prepend(descriptions,prepend)
% [DESCRIPTIONS] = DESCRIPTION_PREPEND(DESCRIPTIONS,PREPEND)
%
%     Prepends each string in the cell array DESCRIPTIONS with the descriptor
%     string PREPEND.

for q = 1:length(descriptions)
  descriptions{q} = strcat(prepend, ' ', descriptions{q});
end
