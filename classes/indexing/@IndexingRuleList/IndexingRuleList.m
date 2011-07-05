classdef IndexingRuleList < Singleton
% IndexingRuleList -- Singleton list of IndexingRule subclasses
%
% self = IndexingRuleList()
%
%     A functionally vacuous class that generates a list of subclasses of
%     IndexingRule. This is done once per startup of Matlab, and this class
%     assumes that all IndexingRule subclasses are contained in the same folder
%     as the folder containing the class IndexingRule.
  properties(SetAccess=private)
    rules = {}; % The cell array list of rules
  end
  methods(Static)
    function self = instance()
      persistent obj
      if isempty(obj)
        obj = IndexingRuleList();
      end
      self = obj;
    end
  end
  methods(Access=private)
    function self = IndexingRuleList()
      location = which('IndexingRule');
      location = fileparts(location); % Discard filename
      fileseplocation = find(location==filesep, 1, 'last'); 
      location = location(1:fileseplocation-1); % Discard @IndexingRule

      % Get class list
      temp = what(location);
      temp = temp.classes;

      % Weed out IndexingRule, IndexingRuleList
      tf = strcmp('IndexingRule', temp);
      temp(tf) = [];
      tf = strcmp('IndexingRuleList', temp);
      temp(tf) = [];

      for q = 1:length(temp)
        self.rules{q} = eval([temp{q} '.instance()']);
      end
    end
  end
end
