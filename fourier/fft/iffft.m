function[nodes,ks] = iffft(modes,varargin);
% [nodes] = iffft(modes,{gamma=0,delta=0,shift=0,scale=1})
%
%     A wrapper for Matlab's IFFT. For gamma=delta=0, takes in modal
%     coefficients corresponding to basis functions 1/sqrt(2*pi)*exp(i*k*theta)
%     and returns the nodal evaluations at the canonical Fourier nodal
%     locations. 
%
%     The inputs gamma and delta specify the type of generalized Fourier
%     expansion, and the values shift and scale dictate the affine map from
%     [-pi,pi] to another interval.

global packages;
opt = packages.common.input_schema({'gamma','delta','shift','scale'},{0,0,0,1},[],varargin{:});
conn = packages.speclab.fourier.connection.negative_integer_separation_connection;
N = length(modes);

if (opt.gamma+opt.delta)>0
  modes = conn(modes,opt.gamma,opt.delta,opt.gamma,opt.delta);
end

ks = packages.speclab.common.integer_range(N);
phase = exp(-i*ks*pi/N);
phase(ks==0) = 1;
phase(mod(ks,2)==1) = phase(mod(ks,2)==1)*-1;
modes = modes./phase;

modes = ifftshift(modes)*N/sqrt(2*pi);
nodes = ifft(modes);
