function[obj] = indexing_parser(self,inp)
% indexing_parser -- Returns an IndexingRule subclass
%
% obj = indexing_parser(input)
% 
%     The given input is a string (or an IndexingRule subclass) that specifies
%     some kind of known indexing rule. This function parses the input and
%     returns the appropriate IndexingRule subclass instance, if appropriate. 
%
%     If the input type is not recognized, an error is NOT thrown -- instead, 
%     the class instance's 'default_indexing_rule' property is returned.

persistent rules
if isempty(rules)
  rules = IndexingRuleList.instance();
  rules = rules.rules;
end

if not(isa(inp, 'char'))
  if isa(inp, 'IndexingRule')
    obj = inp;
    return
  elseif isempty(inp)
    obj = self.default_indexing_rule;
    return
  elseif isa(inp, 'numeric')
    inp = num2str(inp);
  else
    error('I don''t know what to do with this non-character indexing rule');
  end
end

obj = [];
bool = false;
q = 1;
while not(bool) & (q <= length(rules))
  curr_norm = rules{q};
  [bool, obj] = curr_norm.id_compare(inp);
  q = q + 1;
end
if not(bool)
  obj = self.default_indexing_rule;
end
