function[I] = indexing_parser(inp)
% idnexing_parser -- Returns IndexingRule-type object based on input
%
% I = indexing_parser(inp)
%
%     Based on the string id `inp', this function identifies a known indexing
%     rule from the id and returns it. If `inp' is a scalar, it is run through
%     Matlab's num2str before the id process is begun.

persistent rule_list
if isempty(rule_list)
  rule_list = {ZeroBasedIndexing.instance(), ...
               OneBasedIndexing.instance(), ...
               IntegerIndexing.instance()};
               %MultinomialIndexing.instance() };
end

for q = 1:length(rule_list)
  I = rule_list{q};
  tf = I.id_compare(inp);
  if tf
    I = rule_list{q};
    break;
  end
end
