function[b] = second_difference_sobolev_proportionality(alpha)
% second_difference_sobolev_proportionality -- constant of proportionality
%
% b = second_difference_sobolev_proportionality(alpha)
%
%     Calculates and returns the constant of proportionality
%     connecting the central second-difference Stein integral and the order
%     alpha Sobolev seminorm.

exponent = (2*alpha + 1);
if abs(exponent-2)<1e-10
  b = pi/4;
elseif abs(exponent-3)<1e-10
  b = log(2);
elseif abs(exponent-4)<1e-10
  b = pi/3;
else
  % All hail Wolfram Alpha
  b = 2^(exponent-5)*(2^exponent-8)*sin((pi*exponent)/2)*gamma(1-exponent);
end
b = 1./(4*b);

% WTF is this factor about?
b = (2)^(2*alpha-2)*b;
%b = (1/pi)^(2*alpha-2)*b;
