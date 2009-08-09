function[f] = chebifft_online(F,data);
% [f] = chebifft(F,data);
%
%     Online computations for the Chebyshev IFFT. 

if strcmpi(data.points, 'gq')

  f = spdiags(shift,0,data.N,data.N)*F;
  f = [0; real(flipud(f(2:end))); real(f)] + i*[0; -imag(flipud(f(2:end))); imag(f)];
  f = ifft(ifftshift(f,1),[],1)*2*data.N;
  f = real(f(1:data.N));

elseif strcmpi(data.points, 'grq')
  error('Not yet implemented')
elseif strcmpi(data.points, 'glq')
  error('Not yet implemented')
else
  error('Unrecognized flag for location of nodes')
end
