function validate(self, inp)
% validate -- Input validate method for WholeMultiIndices
%
% validate(self, inp)
%
%     Checks if input is an array of multi-indices. (Type checking not
%     performed.) 

assert(size(inp,1)==self.dim, ['Input does not have ' num2str(self.dim) ' rows']);
assert(all(inp(:)>=0), 'Input array does not contain whole numbers')
