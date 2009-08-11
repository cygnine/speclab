function[container] = fourier_connection_tests(container,opt)
% [container] = fourier_connection_tests(container,opt)
%
%     Appends function connection tests for Fourier Series to the
%     TestContainer container. The Fourier family to use is specified by the
%     struct of parameters opt.

import debug.*

test = ValidationTest('description', 'Fourier Positive connection',...
                      'parameters', opt,...
                      'validator', @pconnection_validator,...
                      'data_generator', @pconnection_data);
container = container.append(test);

test = ValidationTest('description', 'Fourier Negative connection',...
                      'parameters', opt,...
                      'validator', @nconnection_validator,...
                      'data_generator', @pconnection_data);
container = container.append(test);

test = ValidationTest('description', 'Fourier online Positive connection',...
                      'parameters', opt,...
                      'validator', @pconnection_online_validator,...
                      'data_generator', @pconnection_online_data);
container = container.append(test);

test = ValidationTest('description', 'Fourier online Negative connection',...
                      'parameters', opt,...
                      'validator', @nconnection_online_validator,...
                      'data_generator', @pconnection_online_data);
container = container.append(test);

end

function[data] = pconnection_data(opt)
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
  v = fourier.eval.fseries(x,ks,opt);
  demoted_modes = v'*spdiags(w,0,opt.N,opt.N)*f(x);

  [data.modes,data.demoted_modes] = deal(modes,demoted_modes);
end

function[tf] = pconnection_validator(data,opt)
  global handles;
  fourier = handles.speclab.fourier;
  pconnect = fourier.connection.positive_integer_separation_connection;

  [modes,demoted_modes] = deal(data.modes,data.demoted_modes);

  modes_promoted = pconnect(demoted_modes,0,0,opt.gamma,opt.delta);

  tol = 10^(-8+(opt.gamma+opt.delta)/4);
  GD = opt.gamma+opt.delta;

  tf = norm(modes_promoted(1+GD:end-GD) - modes(1+GD:end-GD))<tol;
end

function[tf] = nconnection_validator(data,opt)
  global handles;
  fourier = handles.speclab.fourier;
  nconnect = fourier.connection.negative_integer_separation_connection;

  [modes,demoted_modes] = deal(data.modes,data.demoted_modes);

  modes_connect = nconnect(modes,opt.gamma,opt.delta,opt.gamma,opt.delta);

  tol = 10^(-8+(opt.gamma+opt.delta)/4);
  GD = opt.gamma+opt.delta;

  tf = norm(demoted_modes(1+GD:end-GD) - modes_connect(1+GD:end-GD))<tol;
end

function[data] = pconnection_online_data(opt)
  global handles;
  fourier = handles.speclab.fourier;
  ks = handles.speclab.common.integer_range(opt.N);
  sss = handles.speclab.common.standard_scaleshift_1d;

  f = @(x) exp(sin(sss(x,opt)));

  [x,w] = fourier.quad.gauss_quadrature(opt.N,opt);
  v = fourier.eval.fseries(x,ks,opt);
  modes = v'*spdiags(w,0,opt.N,opt.N)*f(x);

  conndata = fourier.connection.integer_separation_connection_overhead(... 
      opt.N,0,0,opt.gamma,opt.delta);

  opt.gamma = 0; opt.delta = 0;
  [x,w] = fourier.quad.gauss_quadrature(opt.N,opt);
  v = fourier.eval.fseries(x,ks,opt);
  demoted_modes = v'*spdiags(w,0,opt.N,opt.N)*f(x);

  [data.modes,data.demoted_modes, data.conndata] = ...
      deal(modes,demoted_modes,conndata);
end

function[tf] = pconnection_online_validator(data,opt)
  global handles;
  fourier = handles.speclab.fourier;
  pconnect = fourier.connection.positive_integer_separation_connection_online;

  [modes,demoted_modes,conndata] =...
    deal(data.modes,data.demoted_modes,data.conndata);

  modes_promoted = pconnect(demoted_modes,conndata);

  tol = 10^(-8+(opt.gamma+opt.delta)/4);
  GD = opt.gamma+opt.delta;

  tf = norm(modes_promoted(1+GD:end-GD) - modes(1+GD:end-GD))<tol;
end

function[tf] = nconnection_online_validator(data,opt)
  global handles;
  fourier = handles.speclab.fourier;
  nconnect = fourier.connection.negative_integer_separation_connection_online;

  [modes,demoted_modes,conndata] = ...
      deal(data.modes,data.demoted_modes,data.conndata);

  modes_connect = nconnect(modes,conndata);

  tol = 10^(-8+(opt.gamma+opt.delta)/4);
  GD = opt.gamma+opt.delta;

  tf = norm(demoted_modes(1+GD:end-GD) - modes_connect(1+GD:end-GD))<tol;
end
