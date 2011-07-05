function validate(self, inp)
% validate -- Input validate method for IntegerNumbers
%
% validate(self, inp)
%
%     Checks if input is an array of integers. (Type checking not performed.)

assert(isa(inp, 'numeric'), 'Input is not an array of integers');
