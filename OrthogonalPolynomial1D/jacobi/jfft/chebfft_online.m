function[F] = chebfft_online(f,data);
% [F] = chebfft_online(f,data)
%
%     Computes the Chebyshev FFT of f using data, the output of
%     chebfft_overhead.

if strcmpi(data.points, 'gq')
  [N,shift] = deal(data.N,data.shift);
  ft = [f; flipud(f)];

  ft = fftshift(fft(ft,[],1),1);

  ft = ft((N+1):(2*N),:);
  
  ft = spdiags(shift,0,N,N)*ft;

  F = real(ft);
elseif strcmpi(opt.points, 'grq')
  error('Not yet implemented')
elseif strcmpi(opt.points, 'glq')
  error('Not yet implemented')
else
  error('Unrecognized flag for location of nodes')
end
