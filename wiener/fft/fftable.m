function[tf,varargout] = fftable(varargin)
% [tf,{S,T}] = fftable({s=1, t=0})
% 
%     This function is meant to take in a struct of parameters specifying a
%     Wiener function family and output a boolean describing whether or not
%     this family can use the FFT. Optionally, it returns the integers A and B
%     denoting the parameter separation of alpha and beta from 1 and 0
%     respectively. 
%
%     A Wiener expansion can support the FFT if s>1/2, t>-1/2, and both s and t are
%     integers.

global handles;
opt = handles.common.InputSchema({'s', 't'}, {1, 0}, [], varargin{:});

tol = 1e-12;
S = opt.s - 1;
T = opt.t;

tf = (abs(S-round(S))<tol) && (abs(T - round(T))<tol);
if tf
  S = round(S);
  T = round(T);
else
  S = [];
  T = [];
end

varargout{1} = S;
varargout{2} = T;
