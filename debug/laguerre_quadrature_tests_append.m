function[container] = laguerre_quadrature_tests_append(container,opt)
% [container] = laguerre_quadrature_tests_append(container,opt)
%
%     Called by laguerre_tests to add ValidationTests given a base laguerre
%     polynomial family descriptor OPT.

import debug.*

test = ValidationTest('description', 'Gauss quadrature nodes inside interval',...
                      'parameters', opt,...
                      'validator', @gauss_quad_nodes,...
                      'data_generator', @gauss_quad_data);
container = container.append(test);

test = ValidationTest('description', 'Gauss quadrature weights positive',...
                      'parameters', opt,...
                      'validator', @gauss_quad_weights,...
                      'data_generator', @gauss_quad_data);
container = container.append(test);

test = ValidationTest('description', 'Gauss quadrature integration accuracy',...
                      'parameters', opt,...
                      'validator', @gauss_quad_acc,...
                      'data_generator', @gauss_quad_acc_data);
container = container.append(test);

opt.r = opt.shift;
test = ValidationTest('description', 'Gauss-Radau quadrature nodes inside interval',...
                      'parameters', opt,...
                      'validator', @gauss_quad_nodes,...
                      'data_generator', @gauss_radau_data);
container = container.append(test);

test = ValidationTest('description', 'Gauss-Radau quadrature weights positive',...
                      'parameters', opt,...
                      'validator', @gauss_quad_weights,...
                      'data_generator', @gauss_radau_data);
container = container.append(test);

test = ValidationTest('description', 'Gauss-Radau quadrature integration accuracy',...
                      'parameters', opt,...
                      'validator', @gauss_quad_acc,...
                      'data_generator', @gauss_radau_quad_acc_data);
container = container.append(test);

function[data] = gauss_quad_data(opt)
  global handles;
  lag = handles.speclab.orthopoly1d.laguerre;

  [x,w] = lag.quad.gauss_quadrature(opt.N,opt);
  data.x = x;
  data.w = w;

function[data] = gauss_radau_data(opt)
  global handles;
  lag = handles.speclab.orthopoly1d.laguerre;

  [x,w] = lag.quad.gauss_radau_quadrature(opt.N,opt);
  data.x = x;
  data.w = w;

function[tf] = gauss_quad_nodes(data,opt)
  global handles;
  lag = handles.speclab.orthopoly1d.laguerre;

  lint = lag.interval(opt);
  tol = 1e-12;
  [x,w] = deal(data.x,data.w);
  tf = all(x<lint(2)+tol) & all(x>lint(1)-tol);

function[tf] = gauss_quad_weights(data,opt)
  global handles;
  lag = handles.speclab.orthopoly1d.laguerre;

  tol = 1e-12;
  [x,w] = deal(data.x,data.w);
  tf = all(w>-tol);

function[data] = gauss_quad_acc_data(opt)
  global handles;
  lag = handles.speclab.orthopoly1d.laguerre;

  [x,w] = lag.quad.gauss_quadrature(opt.N,opt);
  ps = lag.eval.eval_laguerre_poly(x,0:(2*opt.N-1),opt);
  [data.x,data.w,data.ps] = deal(x,w,ps);

function[data] = gauss_radau_quad_acc_data(opt)
  global handles;
  lag = handles.speclab.orthopoly1d.laguerre;

  [x,w] = lag.quad.gauss_radau_quadrature(opt.N,opt);
  ps = lag.eval.eval_laguerre_poly(x,0:(2*opt.N-2),opt);
  [data.x,data.w,data.ps] = deal(x,w,ps);

function[tf] = gauss_quad_acc(data,opt)
  tol = 1e-8;
  [x,w,ps] = deal(data.x,data.w,data.ps);

  intval = zeros(size(ps,2));
  for n = 2:(size(ps,2))
    intval(n) = w'*ps(:,n);
  end

  tf = all(abs(intval)<tol);
