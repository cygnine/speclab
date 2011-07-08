classdef NaturalMultiIndices < IndexSet
% NaturalMultiIndices -- The set of natural-number-valued multi-indices of a given dimension
%
% self = NaturalMultiIndices(dim)
%
%     The set of multi-indices in dim-dimensional Euclidean space with
%     values in NaturalNumbers for each dimension. This is a class with a
%     singleton instance for each dim.
%
%     An element of this set is a column vector with dim entries. An array
%     is a collection of these elements if it has dim rows.

  properties(SetAccess=private)
    descriptor
    dim
  end
  methods(Static)
    function self = instance(dim)
      persistent objs
      if nargin < 1
        dim = 1;
      end
      if isempty(objs)
        objs = {};
      end
      if (length(objs) < dim) || isempty(objs{dim})
        objs{dim} = NaturalMultiIndices(dim);
      end
      self = objs{dim};
    end
  end
  methods(Access=private)
    function self = NaturalMultiIndices(dim)
      self.dim = dim;
      self.descriptor = ['Natural-number-valued, ' num2str(dim) '-dimensional multi-indices'];
    end
  end
  methods
    validate(self,inp)
  end
end
