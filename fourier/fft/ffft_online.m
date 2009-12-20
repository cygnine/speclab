function[modes] = ffft_online(nodes,data);
% [modes] = ffft_online(nodes,data)
% 
%     The online "Fourier FFT".

persistent connection
if isempty(connection)
  from speclab.fourier.connection import positive_integer_separation_connection_online as connection
end

modes = fftshift(fft(nodes));

modes = data.phase.*modes;

if data.GD>0
  modes = connection(modes,data.conndata);
end
