function[tf,A,B] = fftable(varargin)
% [tf,A,B] = fftable({alpha=-1/2, beta=-1/2})
% 
%     This function is meant to take in a struct of parameters specifying a
%     Jacobi polynomial family and output a boolean describing whethet or not
%     this family can use the FFT. Optionally, it returns the integers A and B
%     denoting the parameter separation of alpha and beta from -1/2,
%     respectively. 

global packages;
opt = packages.common.input_schema({'alpha', 'beta'}, {-1/2, -1/2}, [], varargin{:});

tol = 1e-12;
A = opt.alpha + 1/2;
B = opt.beta + 1/2;

tf = (abs(A-round(A))<tol) && (abs(B - round(B))<tol);
if tf
  A = round(A);
  B = round(B);
else
  A = [];
  B = [];
end
