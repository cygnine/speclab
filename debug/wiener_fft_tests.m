function[container] = wiener_fft_tests(container,opt)
% [container] = wiener_fft_tests(container,opt)
%
%     Appends function fft tests for Wiener functions to the TestContainer
%     container. The Wiener family to use is specified by the struct of
%     parameters opt.

import debug.*

test = ValidationTest('description', 'Wiener FFT collocation',...
                      'parameters', opt,...
                      'validator', @fftcoll_validator,...
                      'data_generator', @fftcoll_data);
container = container.append(test);

test = ValidationTest('description', 'Wiener IFFT collocation',...
                      'parameters', opt,...
                      'validator', @ifftcoll_validator,...
                      'data_generator', @ifftcoll_data);
container = container.append(test);

test = ValidationTest('description', 'online Wiener FFT collocation',...
                      'parameters', opt,...
                      'validator', @online_fftcoll_validator,...
                      'data_generator', @online_fftcoll_data);
container = container.append(test);

test = ValidationTest('description', 'online Wiener IFFT collocation',...
                      'parameters', opt,...
                      'validator', @online_ifftcoll_validator,...
                      'data_generator', @online_ifftcoll_data);
container = container.append(test);

test = ValidationTest('description', 'Wiener FFT Galerkin',...
                      'parameters', opt,...
                      'validator', @fftgal_validator,...
                      'data_generator', @fftgal_data);
container = container.append(test);

test = ValidationTest('description', 'Wiener IFFT Galerkin',...
                      'parameters', opt,...
                      'validator', @ifftgal_validator,...
                      'data_generator', @ifftgal_data);
container = container.append(test);

test = ValidationTest('description', 'online Wiener FFT Galerkin',...
                      'parameters', opt,...
                      'validator', @online_fftgal_validator,...
                      'data_generator', @online_fftgal_data);
container = container.append(test);

test = ValidationTest('description', 'online Wiener IFFT Galerkin',...
                      'parameters', opt,...
                      'validator', @online_ifftgal_validator,...
                      'data_generator', @online_ifftgal_data);
container = container.append(test);

function[data] = fftcoll_data(opt)
  global packages;
  wiener = packages.speclab.wiener;

  sss = packages.speclab.common.standard_scaleshift_1d;
  f = @(x) exp(-sss(x,opt).^2)./(1+sss(x,opt).^2);
  sopt = opt; sopt.s = 1; sopt.t = 0;

  [x,w] = wiener.quad.pi_gauss_quadrature(opt.N,sopt);
  ks = packages.speclab.common.integer_range(opt.N);
  ws = wiener.eval.wiener_function(x,ks,opt);
  fx = f(x);
  modes = ws'*(fx.*w);

  fft_modes = wiener.fft.wfft_collocation(fx,opt);

  [data.modes, data.fft_modes] = deal(modes, fft_modes);

function[tf] = fftcoll_validator(data,opt)

  [modes,fft_modes] = deal(data.modes, data.fft_modes);

  tol = 10^(-6+opt.s/5);
  ST = opt.s + opt.t + 1;

  tf = norm(modes(ST:(end-ST-1)) - fft_modes(ST:(end-ST-1)))<tol;

function[data] = ifftcoll_data(opt)
  global packages;
  wiener = packages.speclab.wiener;

  sss = packages.speclab.common.standard_scaleshift_1d;
  f = @(x) exp(-sss(x,opt).^2)./(1+sss(x,opt).^2);
  sopt = opt; sopt.s = 1; sopt.t = 0;

  [x,w] = wiener.quad.pi_gauss_quadrature(opt.N,sopt);
  ks = packages.speclab.common.integer_range(opt.N);
  ws = wiener.eval.wiener_function(x,ks,opt);
  fx = f(x);
  modes = ws'*(fx.*w);

  fft_modes = wiener.fft.wfft_collocation(fx,opt);

  [data.modes, data.fx] = deal(modes, fx);

function[tf] = ifftcoll_validator(data,opt)
  
  global packages;
  wiener = packages.speclab.wiener;
  
  [modes,fx] = deal(data.modes, data.fx);

  fx_fft = wiener.fft.iwfft_collocation(modes,opt);

  tol = 10^(-6+opt.s/5);

  tf = norm(fx_fft - fx)<tol;

function[data] = online_fftcoll_data(opt)
  global packages;
  wiener = packages.speclab.wiener;

  sss = packages.speclab.common.standard_scaleshift_1d;
  f = @(x) exp(-sss(x,opt).^2)./(1+sss(x,opt).^2);
  sopt = opt; sopt.s = 1; sopt.t = 0;

  [x,w] = wiener.quad.pi_gauss_quadrature(opt.N,sopt);
  fx = f(x);
  fft_modes = wiener.fft.wfft_collocation(fx,opt);

  fft_data = wiener.fft.wfft_collocation_overhead(opt.N,opt);

  [data.fft_modes, data.fx, data.fft_data] = deal(fft_modes, fx, fft_data);

function[tf] = online_fftcoll_validator(data,opt)
  global packages;
  wiener = packages.speclab.wiener;

  [fft_modes,fx,fft_data] = deal(data.fft_modes, data.fx, data.fft_data);
  fft_modes_online = wiener.fft.wfft_collocation_online(fx,fft_data);

  tol = 1e-6;
  tf = norm(fft_modes - fft_modes_online)<tol;

function[data] = online_ifftcoll_data(opt)
  global packages;
  wiener = packages.speclab.wiener;

  sss = packages.speclab.common.standard_scaleshift_1d;
  f = @(x) exp(-sss(x,opt).^2)./(1+sss(x,opt).^2);
  sopt = opt; sopt.s = 1; sopt.t = 0;

  [x,w] = wiener.quad.pi_gauss_quadrature(opt.N,sopt);
  fx = f(x);
  fft_modes = wiener.fft.wfft_collocation(fx,opt);
  fft_data = wiener.fft.wfft_collocation_overhead(opt.N,opt);

  [data.fft_modes, data.fx, data.fft_data] = ...
    deal(fft_modes, fx, fft_data);

function[tf] = online_ifftcoll_validator(data,opt)
  global packages;
  wiener = packages.speclab.wiener;

  [fft_modes, fx, fft_data] = deal(data.fft_modes, data.fx, data.fft_data);
  fx_online = wiener.fft.iwfft_collocation_online(fft_modes,fft_data);

  tol = 1e-6;
  tf = norm(fx-fx_online)<tol;

function[data] = fftgal_data(opt)
  global packages;
  wiener = packages.speclab.wiener;

  sss = packages.speclab.common.standard_scaleshift_1d;
  f = @(x) exp(-sss(x,opt).^2)./(1+sss(x,opt).^2);
  sopt = opt; sopt.s = 1; sopt.t = 0;

  [x,w] = wiener.quad.pi_gauss_quadrature(opt.N,sopt);
  ks = packages.speclab.common.integer_range(opt.N);
  ws = wiener.eval.wiener_function(x,ks,opt);
  fx = f(x);
  modes = ws'*(fx.*w);

  fft_modes = wiener.fft.wfft_galerkin(fx,opt);

  [data.modes, data.fft_modes] = deal(modes, fft_modes);

function[tf] = fftgal_validator(data,opt)

  [modes,fft_modes] = deal(data.modes, data.fft_modes);

  tol = 10^(-6+opt.s/5);
  ST = opt.s + opt.t + 1;

  tf = norm(modes(ST:(end-ST-1)) - fft_modes(ST:(end-ST-1)))<tol;

function[data] = ifftgal_data(opt)
  global packages;
  wiener = packages.speclab.wiener;

  sss = packages.speclab.common.standard_scaleshift_1d;
  f = @(x) exp(-sss(x,opt).^2)./(1+sss(x,opt).^2);
  sopt = opt; sopt.s = 1; sopt.t = 0;

  [x,w] = wiener.quad.pi_gauss_quadrature(opt.N,sopt);
  ks = packages.speclab.common.integer_range(opt.N);
  ws = wiener.eval.wiener_function(x,ks,opt);
  fx = f(x);
  modes = ws'*(fx.*w);

  [data.modes, data.fx] = deal(modes, fx);

function[tf] = ifftgal_validator(data,opt)
  
  global packages;
  wiener = packages.speclab.wiener;
  
  [modes,fx] = deal(data.modes, data.fx);

  fx_fft = wiener.fft.iwfft_galerkin(modes,opt);

  tol = 10^(-6+opt.s/5);

  tf = norm(fx_fft - fx)<tol;

function[data] = online_fftgal_data(opt)
  global packages;
  wiener = packages.speclab.wiener;

  sss = packages.speclab.common.standard_scaleshift_1d;
  f = @(x) exp(-sss(x,opt).^2)./(1+sss(x,opt).^2);
  sopt = opt; sopt.s = 1; sopt.t = 0;

  [x,w] = wiener.quad.pi_gauss_quadrature(opt.N,sopt);
  fx = f(x);
  fft_modes = wiener.fft.wfft_galerkin(fx,opt);

  fft_data = wiener.fft.wfft_galerkin_overhead(opt.N,opt);

  [data.fft_modes, data.fx, data.fft_data] = deal(fft_modes, fx, fft_data);

function[tf] = online_fftgal_validator(data,opt)
  global packages;
  wiener = packages.speclab.wiener;

  [fft_modes,fx,fft_data] = deal(data.fft_modes, data.fx, data.fft_data);
  fft_modes_online = wiener.fft.wfft_galerkin_online(fx,fft_data);

  tol = 1e-6;
  tf = norm(fft_modes - fft_modes_online)<tol;

function[data] = online_ifftgal_data(opt)
  global packages;
  wiener = packages.speclab.wiener;

  sss = packages.speclab.common.standard_scaleshift_1d;
  f = @(x) exp(-sss(x,opt).^2)./(1+sss(x,opt).^2);
  sopt = opt; sopt.s = 1; sopt.t = 0;

  [x,w] = wiener.quad.pi_gauss_quadrature(opt.N,sopt);
  fx = f(x);
  fft_modes = wiener.fft.wfft_galerkin(fx,opt);
  fft_data = wiener.fft.wfft_galerkin_overhead(opt.N,opt);

  [data.fft_modes, data.fx, data.fft_data] = ...
    deal(fft_modes, fx, fft_data);

function[tf] = online_ifftgal_validator(data,opt)
  global packages;
  wiener = packages.speclab.wiener;

  [fft_modes, fx, fft_data] = deal(data.fft_modes, data.fx, data.fft_data);
  fx_online = wiener.fft.iwfft_galerkin_online(fft_modes,fft_data);

  tol = 1e-6;
  tf = norm(fx-fx_online)<tol;
