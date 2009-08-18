function[container] = fourier_fft_tests(container,opt)
% [container] = fourier_fft_tests(container,opt)
%
%     Appends fft tests for Fourier Series to the TestContainer container. The
%     Fourier family to use is specified by the struct of parameters opt.

import debug.*

test = ValidationTest('description', 'Fourier fft',...
                      'parameters', opt,...
                      'validator', @fft_validator,...
                      'data_generator', @fft_data);
container = container.append(test);

test = ValidationTest('description', 'Fourier ifft',...
                      'parameters', opt,...
                      'validator', @ifft_validator,...
                      'data_generator', @fft_data);
container = container.append(test);

test = ValidationTest('description', 'Fourier fft online',...
                      'parameters', opt,...
                      'validator', @fft_online_validator,...
                      'data_generator', @fft_online_data);
container = container.append(test);

test = ValidationTest('description', 'Fourier ifft online',...
                      'parameters', opt,...
                      'validator', @ifft_online_validator,...
                      'data_generator', @fft_online_data);
container = container.append(test);

function[data] = fft_data(opt)
  global handles;
  fourier = handles.speclab.fourier;
  ks = handles.speclab.common.integer_range(opt.N);
  sss = handles.speclab.common.standard_scaleshift_1d;

  f = @(x) exp(sin(sss(x,opt)));

  [x,w] = fourier.quad.gauss_quadrature(opt.N,opt);
  v = fourier.eval.fseries(x,ks,opt);
  modes = v'*spdiags(w,0,opt.N,opt.N)*f(x);

  opt.gamma = 0; opt.delta = 0;
  [x,w] = fourier.quad.gauss_quadrature(opt.N,opt);
  fx = f(x);

  [data.modes,data.fx] = deal(modes,fx);

function[tf] = fft_validator(data,opt)
  global handles;
  fourier = handles.speclab.fourier;

  [modes,fx] = deal(data.modes,data.fx);

  tol = 10^(-8+(opt.gamma+opt.delta)/4);

  fft_modes = fourier.fft.ffft(fx,opt);

  GD = opt.gamma+opt.delta;

  tf = norm(fft_modes(1+GD:end-GD) - modes(1+GD:end-GD))<tol;

function[tf] = ifft_validator(data,opt)
  global handles;
  fourier = handles.speclab.fourier;

  [modes,fx] = deal(data.modes,data.fx);

  tol = 10^(-8+(opt.gamma+opt.delta)/4);

  ifft_fx= fourier.fft.iffft(modes,opt);

  tf = norm(fx - ifft_fx)<tol;

function[data] = fft_online_data(opt)
  global handles;
  fourier = handles.speclab.fourier;
  ks = handles.speclab.common.integer_range(opt.N);
  sss = handles.speclab.common.standard_scaleshift_1d;
  fft_overhead = fourier.fft.ffft_overhead;

  f = @(x) exp(sin(sss(x,opt)));

  [x,w] = fourier.quad.gauss_quadrature(opt.N,opt);
  v = fourier.eval.fseries(x,ks,opt);
  modes = v'*spdiags(w,0,opt.N,opt.N)*f(x);

  fftdata = fft_overhead(opt.N,opt);

  opt.gamma = 0; opt.delta = 0;
  [x,w] = fourier.quad.gauss_quadrature(opt.N,opt);
  fx = f(x);

  [data.modes,data.fx,data.fftdata] = deal(modes,fx,fftdata);

function[tf] = fft_online_validator(data,opt)
  global handles;
  fourier = handles.speclab.fourier;

  [modes,fx,fftdata] = deal(data.modes,data.fx,data.fftdata);

  tol = 10^(-8+(opt.gamma+opt.delta)/4);

  fft_modes = fourier.fft.ffft_online(fx,fftdata);

  GD = opt.gamma+opt.delta;

  tf = norm(fft_modes(1+GD:end-GD) - modes(1+GD:end-GD))<tol;

function[tf] = ifft_online_validator(data,opt)
  global handles;
  fourier = handles.speclab.fourier;

  [modes,fx,fftdata] = deal(data.modes,data.fx,data.fftdata);

  tol = 10^(-8+(opt.gamma+opt.delta)/4);

  ifft_fx= fourier.fft.iffft_online(modes,fftdata);

  tf = norm(fx - ifft_fx)<tol;
