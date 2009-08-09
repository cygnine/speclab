function self = make_vandermonde_inverse(self)
% SELF = MAKE_VANDERMONDE_INVERSE(SELF)
% 
%     Uses the data for the class instance SELF to construct the inverse of the
%     Vandermonde matrix. This method is not valid in this base class. 

% If vandermonde is not created, create it:
if isempty(self.vandermonde)
  self = self.make_vandermonde;
end

size_v = size(self.vandermonde);
if size_v(1)~=size_v(2)
  error('Cannot invert a non-square vandermonde matrix');
else
  self.vandermonde_inverse = inv(self.vandermonde);
end
