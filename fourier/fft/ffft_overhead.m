function[data] = ffft_overhead(N,varargin);
% [data] = ffft_overhead(N,{gamma=0,delta=0,shift=0,scale=1});
%
%     Overhead computation and storage for the forward and inverse Fourier FFT.

global handles;
opt = handles.common.input_schema({'gamma','delta','shift','scale'},{0,0,0,1},[],varargin{:});
conn = handles.speclab.fourier.connection.integer_separation_connection_overhead;

conndata = conn(N,0,0,opt.gamma,opt.delta);

ks = handles.speclab.common.integer_range(N);
phase = exp(-i*ks*pi/N);
phase(ks==0) = 1;
phase(mod(ks,2)==1) = phase(mod(ks,2)==1)*-1;
phase = phase/N*sqrt(2*pi);

GD = opt.gamma+opt.delta;

[data.N, data.phase, data.conndata,data.GD] = deal(N,phase,conndata,GD);
