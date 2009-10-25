function[modes] = ffft_online(nodes,data);
% [modes] = ffft_online(nodes,data)
% 
%     The online "Fourier FFT".

global packages;
conn = packages.speclab.fourier.connection.positive_integer_separation_connection_online;

modes = fftshift(fft(nodes));

modes = data.phase.*modes;

if data.GD>0
  modes = conn(modes,data.conndata);
end
