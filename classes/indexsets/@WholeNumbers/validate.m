function validate(self, inp)
% validate -- Input validate method for WholeNumbers
%
% validate(self, inp)
%
%     Checks if input is an array of whole numbers. (Type checking not
%     performed.)

assert(all(inp(:)>=0), 'Input is not an array of whole numbers')
