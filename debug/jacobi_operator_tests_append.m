function[container] = jacobi_operator_tests_append(container,opt);
% [container] = jacobi_operator_tests_append(container,opt);
%
%     Appends some Jacobi coefficient ValidationTest's to the TestContainer
%     container. The Jacobi family is defined by opt. 

import debug.*
global handles;
jac = handles.speclab.orthopoly1d.jacobi;

test = ValidationTest('description', 'Stiffness operator',...
                      'parameters', opt,...
                      'validator', @stiffness_validator,...
                      'data_generator', @stiffness_data);
container = container.append(test);

[tf,A,B] = jac.jfft.fftable(opt);

if tf
  if (A+B)==0
    test = ValidationTest('description', 'Chebyshev FFT',...
                          'parameters', opt,...
                          'validator', @chebfft_validator,...
                          'data_generator', @chebfft_data);
    container = container.append(test);

    test = ValidationTest('description', 'Chebyshev online FFT',...
                          'parameters', opt,...
                          'validator', @chebfft_online_validator,...
                          'data_generator', @chebfft_online_data);
    container = container.append(test);

    test = ValidationTest('description', 'Chebyshev IFFT',...
                          'parameters', opt,...
                          'validator', @chebifft_validator,...
                          'data_generator', @chebifft_data);
    container = container.append(test);

    test = ValidationTest('description', 'Chebyshev online IFFT',...
                          'parameters', opt,...
                          'validator', @chebifft_online_validator,...
                          'data_generator', @chebifft_online_data);
    container = container.append(test);
  end

  test = ValidationTest('description', 'Jacobi FFT',...
                        'parameters', opt,...
                        'validator', @jfft_validator,...
                        'data_generator', @jfft_data);
  container = container.append(test);

  test = ValidationTest('description', 'Jacobi online FFT',...
                        'parameters', opt,...
                        'validator', @jfft_online_validator,...
                        'data_generator', @jfft_online_data);
  container = container.append(test);

  test = ValidationTest('description', 'Jacobi IFFT',...
                        'parameters', opt,...
                        'validator', @jifft_validator,...
                        'data_generator', @jifft_data);
  container = container.append(test);

  test = ValidationTest('description', 'Jacobi online IFFT',...
                        'parameters', opt,...
                        'validator', @jifft_online_validator,...
                        'data_generator', @jifft_online_data);
  container = container.append(test);
end

function[data] = stiffness_data(opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  modes = randn([opt.N,1]);
  [r,w] = jac.quad.gauss_quadrature(opt.N,opt);
  opt.d = 1;
  dps = jac.eval.eval_jacobi_poly(r,0:(opt.N-1),opt);
  opt.d = 0;
  ps = jac.eval.eval_jacobi_poly(r,0:(opt.N-1),opt);
  vinv = ps'*spdiags(w,0,opt.N,opt.N);

  dfun_rec = dps*modes;
  dmodes = vinv*dfun_rec;

  [data.modes,data.dmodes] = deal(modes,dmodes);

function[tf] = stiffness_validator(data,opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  [modes,dmodes] = deal(data.modes,data.dmodes);

  dmodes2 = jac.operators.stiffness_operator(modes,opt);
  tol = 10^(-7 + opt.alpha/3 + opt.beta/3 + abs(opt.alpha-opt.beta)/3);
  tf = norm(dmodes - dmodes2)<tol;

function[data] = chebfft_data(opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  f = @(x) sin(opt.N/(10*opt.scale)*x);

  [x,w] = jac.quad.gauss_quadrature(opt.N,opt);
  ps = jac.eval.eval_jacobi_poly(x,0:(opt.N-1),opt);
  modes = ps'*spdiags(w,0,opt.N,opt.N)*f(x);
  fx = f(x);

  [data.modes, data.fx] = deal(modes,fx);

function[tf] = chebfft_validator(data,opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  tol = 1e-8;
  [modes,fx] = deal(data.modes,data.fx);
  modes2 = jac.jfft.chebfft(fx,opt);
  tf = norm(modes-modes2)<tol;

function[data] = chebfft_online_data(opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  f = @(x) sin(opt.N/(10*opt.scale)*x);

  [x,w] = jac.quad.gauss_quadrature(opt.N,opt);
  ps = jac.eval.eval_jacobi_poly(x,0:(opt.N-1),opt);
  modes = ps'*spdiags(w,0,opt.N,opt.N)*f(x);
  fx = f(x);
  fftdata = jac.jfft.chebfft_overhead(opt.N,opt);

  [data.modes, data.fx, data.fftdata] = deal(modes,fx,fftdata);

function[tf] = chebfft_online_validator(data,opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  tol = 1e-8;
  [modes,fx,fftdata] = deal(data.modes,data.fx,data.fftdata);
  modes2 = jac.jfft.chebfft_online(fx,fftdata);
  tf = norm(modes-modes2)<tol;

function[data] = chebifft_data(opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  f = @(x) sin(opt.N/(10*opt.scale)*x);
  [x,w] = jac.quad.gauss_quadrature(opt.N,opt);
  modes = jac.jfft.chebfft(f(x),opt);
  [data.fx, data.modes] = deal(f(x),modes);

function[tf] = chebifft_validator(data,opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  tol = 1e-8;
  [fx,modes] = deal(data.fx,data.modes);
  fx2 = jac.jfft.chebifft(modes,opt);
  tf = norm(fx-fx2)<tol;

function[data] = chebifft_online_data(opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  f = @(x) sin(opt.N/(10*opt.scale)*x);
  [x,w] = jac.quad.gauss_quadrature(opt.N,opt);
  modes = jac.jfft.chebfft(f(x),opt);
  fftdata = jac.jfft.chebifft_overhead(opt.N,opt);
  [data.fx, data.modes, data.fftdata] = deal(f(x),modes,fftdata);

function[tf] = chebifft_online_validator(data,opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  tol = 1e-8;
  [fx,modes,fftdata] = deal(data.fx,data.modes,data.fftdata);
  fx2 = jac.jfft.chebifft_online(modes,fftdata);
  tf = norm(fx-fx2)<tol;

function[data] = jfft_data(opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  f = @(x) sin(opt.N/(10*opt.scale)*x);
  [x,w] = jac.quad.gauss_quadrature(opt.N,opt);
  ps = jac.eval.eval_jacobi_poly(x,0:(opt.N-1), opt);
  modes = ps'*spdiags(w,0,opt.N,opt.N)*f(x);

  chebopt = opt;
  chebopt.alpha = -1/2;
  chebopt.beta = -1/2;
  [x,w] = jac.quad.gauss_quadrature(chebopt.N,chebopt);
  fx = f(x);
  [data.modes,data.fx] = deal(modes,fx);

function[tf] = jfft_validator(data,opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  [tf,A,B] = jac.jfft.fftable(opt);
  tol = 10^(-8+(A+B)/4);
  [modes,fx] = deal(data.modes,data.fx);
  modes2 = jac.jfft.jfft(fx,opt);

  tf = norm(modes(1:(end-A-B)) - modes2(1:(end-A-B)))<tol;

function[data] = jfft_online_data(opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  f = @(x) sin(opt.N/(10*opt.scale)*x);
  [x,w] = jac.quad.gauss_quadrature(opt.N,opt);
  ps = jac.eval.eval_jacobi_poly(x,0:(opt.N-1), opt);
  modes = ps'*spdiags(w,0,opt.N,opt.N)*f(x);

  chebopt = opt;
  chebopt.alpha = -1/2;
  chebopt.beta = -1/2;
  [x,w] = jac.quad.gauss_quadrature(chebopt.N,chebopt);
  fx = f(x);
  
  fftdata = jac.jfft.jfft_overhead(opt.N,opt);
  [data.modes,data.fx,data.fftdata] = deal(modes,fx,fftdata);

function[tf] = jfft_online_validator(data,opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  [tf,A,B] = jac.jfft.fftable(opt);
  tol = 10^(-8+(A+B)/4);
  [modes,fx,fftdata] = deal(data.modes,data.fx,data.fftdata);
  modes2 = jac.jfft.jfft_online(fx,fftdata);

  tf = norm(modes(1:(end-A-B)) - modes2(1:(end-A-B)))<tol;

function[data] = jifft_data(opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  f = @(x) sin(opt.N/(10*opt.scale)*x);
  chebopt = opt;
  chebopt.alpha = -1/2;
  chebopt.beta = -1/2;
  [x,w] = jac.quad.gauss_quadrature(chebopt.N,chebopt);

  fx = f(x);
  modes = jac.jfft.jfft(fx,opt);

  [data.fx,data.modes] = deal(fx,modes);

function[tf] = jifft_validator(data,opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  [tf,A,B] = jac.jfft.fftable(opt);
  tol = 10^(-8+(A+B)/4);
  [fx,modes] = deal(data.fx,data.modes);
  fx2 = jac.jfft.jifft(modes,opt);

  tf = norm(fx-fx2)<tol;

function[data] = jifft_online_data(opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  f = @(x) sin(opt.N/(10*opt.scale)*x);
  chebopt = opt;
  chebopt.alpha = -1/2;
  chebopt.beta = -1/2;
  [x,w] = jac.quad.gauss_quadrature(chebopt.N,chebopt);

  fx = f(x);
  modes = jac.jfft.jfft(fx,opt);
  fftdata = jac.jfft.jifft_overhead(opt.N,opt);

  [data.fx,data.modes,data.fftdata] = deal(fx,modes,fftdata);

function[tf] = jifft_online_validator(data,opt);
  global handles;
  jac = handles.speclab.orthopoly1d.jacobi;

  [tf,A,B] = jac.jfft.fftable(opt);
  tol = 10^(-8+(A+B)/4);
  [fx,modes,fftdata] = deal(data.fx,data.modes,data.fftdata);
  fx2 = jac.jfft.jifft_online(modes,fftdata);

  tf = norm(fx-fx2)<tol;
