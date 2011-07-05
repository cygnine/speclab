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
               IntegerIndexing.instance(), ...
               TotalGLexIndexing.instance(), ...
               TotalGRevLexIndexing.instance(), ...
               MarginalGLexIndexing.instance(), ...
               MarginalGRevLexIndexing.instance()};
               %MultinomialIndexing.instance() };
end

% If we're given an IndexingRule instance, just return it
if isa(inp, 'IndexingRule')
  I = inp;
  return
end

% Otherwise search string id's
for q = 1:length(rule_list)
  I = rule_list{q};
  tf = I.id_compare(inp);
  if tf
    break;
  end
  I = [];
end

if isempty(I)
  error('String input is not an ID for an indexing rule');
end
