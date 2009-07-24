function[nodes] = ifft(modes);
% [NODES] = IFFT(MODES,{GAMMA=0,DELTA=0,SHIFT=0,SCALE=1})
%
%     A wrapper for Matlab's IFFT. For GAMMA=DELTA=0, takes in modal
%     coefficients corresponding to basis functions 1/sqrt(2*pi)*exp(i*k*theta)
%     and returns the nodal evaluations at the canonical Fourier nodal
%     locations. 
%
%     The inputs GAMMA and DELTA specify the type of generalized Fourier
%     expansion, and the values SHIFT and SCALE dictate the affine map from
%     [-pi,pi] to another interval.

global handles;
opt = handles.common.InputSchema({'gamma','delta','shift','scale'},{0,0,0,1},[],varargin{:});
conn = handles.speclab.fourier.connection.negative_integer_separation_connection;

if (opt.gamma+opt.delta)>0
  modes = conn(modes,0,0,opt.gamma,opt.delta)
end

modes = ifftshift(modes)*length(modes)/sqrt(2*pi);
nodes = ifft(modes);
