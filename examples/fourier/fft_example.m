%%% Example script for Fourier approximation
%%% FFT: calling syntax for FFT

% Matlab complains bitterly if you create a file fft.m and it sees it in the
% path. To avoid such problems, I've renamed the "fast Fourier transform" as the
% "Fourier fast Fourier transform". ffft instead of fft. iffft instead of ifft.

clear
global packages;
fourier = packages.speclab.fourier;

% The regular FFT is possible:
N = 277;
[theta,w] = fourier.quad.gauss_quadrature(N);
f = @(x) atan(5*sin(6*x));
ft = f(theta);
modes = fourier.fft.ffft(ft);

% Likewise, there's an offline/online decomposition
fftdata = fourier.fft.ffft_overhead(N); % no optional inputs here
modes2 = fourier.fft.ffft_online(ft,fftdata);

% You can compare modes and modes2.

% The FFT may also be used for the gamma, delta integers. I've only tested it
% for gamma different from 0. 
fopt.gamma = 5;
modes3 = fourier.fft.ffft(ft, fopt);
fftdata = fourier.fft.ffft_overhead(N,fopt); % optional inputs since gamma!=0
modes4 = fourier.fft.ffft_online(ft,fftdata);
% When doing an online/overhead computation, the necessary optional inputs are
% already stored in fftdata, so you don't need to input them again.

% Again, you can compare modes3 and modes4

% The iffft's work in exactly the same way, *except* that the overhead
% computation is the same function ffft_overhead. There is no iffft_overhead
% function because the data is the same. So we can use the same fftdata computed
% above to calculate the inverse tranform:
ft2 = fourier.fft.iffft_online(modes4, fftdata);

% And you can check that ft2 and ft are the same.
