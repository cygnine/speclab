function[F] = chebfft_online(f,data);
% [F] = chebfft_online(f,data)
%
%     Computes the Chebyshev FFT of f using data, the output of
%     chebfft_overhead.

if strcmpi(data.points, 'gq')
  [N,shift] = deal(data.N,data.shift);
  %%%%%%%%%%% The FFT way -- almost always slower
  %data.ft = [f; flipud(f)];

  %data.ft = fftshift(fft(data.ft,[],1),1);

  %F = data.ft((N+1):(2*N),:);
  
  %F = spdiags(shift,0,N,N)*F;

  %F = real(F);
  %%%%%%%%%%%

  %%%%%%%%%%% The DCT way -- usually faster
  F = dct(f).*shift;
  %%%%%%%%%%%
elseif strcmpi(opt.points, 'grq')
  error('Not yet implemented')
elseif strcmpi(opt.points, 'glq')
  error('Not yet implemented')
else
  error('Unrecognized flag for location of nodes')
end
