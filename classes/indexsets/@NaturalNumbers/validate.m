function validate(self, inp)
% validate -- Input validate method for NaturalNumbers
%
% validate(self, inp)
%
%     Checks if input is an array of natural numbers. (Type checking not
%     performed.)

assert(all(inp(:)>=1), 'Input is not an array of natural numbers')
