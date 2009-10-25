function[nodes] = iffft_online(modes,data);
% [nodes] = iffft_online(modes,data)
% 
%     The online "inverse Fourier FFT".

global packages;
conn = packages.speclab.fourier.connection.negative_integer_separation_connection_online;

if data.GD>0
  modes = conn(modes,data.conndata);
end

modes = modes./data.phase;

nodes = ifft(ifftshift(modes));
