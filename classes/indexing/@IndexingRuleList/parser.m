function[rule] = parser(self, id, default)
% parser -- Returns an IndexingRule subclass
%
% obj = parser(self, id)
% 
%     The given input 'id' is a string (or an IndexingRule subclass) that specifies
%     some kind of known indexing rule. This function parses the input and
%     returns the appropriate IndexingRule subclass instance, if appropriate.
%     This is done by comparing the given id using the compare_id method of each
%     rule in self.rules. If no matching rule is found, an error is thrown.
%
% obj = parser(self, id, default)
%
%     If a IndexingRule subclass instance 'default' is given, it is returned
%     when no matching id in self.rules is found.

if not(isa(id, 'char'))
  if isa(id, 'IndexingRule')
    rule = id;
    return
  elseif isa(id, 'numeric')
    id = num2str(id);
  else
    error('I don''t know what to do with this non-character, non-IndexingRule id');
  end
end

rule = [];
bool = false;
q = 1;
while not(bool) && (q <= length(self.rules))
  curr_norm = self.rules{q};
  [bool, rule] = curr_norm.id_compare(id);
  q = q + 1;
end
if not(bool)
  if nargin < 3
    error('Given id does not match any known indexing rules');
  elseif not(isa(default, 'IndexingRule'));
    error('Given default indexing rule is not an IndexingRule instance');
  else
    rule = default;
  end
end
