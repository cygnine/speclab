function[nodes] = iffft_online(modes,data);
% [nodes] = iffft_online(modes,data)
% 
%     The online "inverse Fourier FFT".

persistent connection
if isempty(connection)
  from speclab.fourier.connection import negative_integer_separation_connection_online as connection
end

if data.GD>0
  modes = connection(modes,data.conndata);
end

modes = modes./data.phase;

nodes = ifft(ifftshift(modes));
