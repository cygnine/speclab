function[container] = wiener_approximation_tests(container,opt)
% [container] = wiener_approximation_tests(container,opt)
%
%     Appends function approximation tests for Wiener functions to the
%     TestContainer container. The Wiener family to use is specified by the
%     struct of parameters opt.

import debug.*

test = ValidationTest('description', 'Mass matrix',...
                      'parameters', opt,...
                      'validator', @mass_validator,...
                      'data_generator', @mass_data);
container = container.append(test);

test = ValidationTest('description', 'Wiener interpolation',...
                      'parameters', opt,...
                      'validator', @interpolation_validator,...
                      'data_generator', @interpolation_data);
container = container.append(test);

test = ValidationTest('description', 'Wiener derivative',...
                      'parameters', opt,...
                      'validator', @derivative_validator,...
                      'data_generator', @derivative_data);
container = container.append(test);

function[data] = mass_data(opt)
  
  from speclab import wiener
  from speclab.common import integer_range
  %wiener = packages.speclab.wiener;
  [x,w] = wiener.quad.pi_gauss_quadrature(opt.N,opt);
  ks = integer_range(opt.N);
  ws = wiener.eval.wiener_function(x,ks,opt);

  mass = ws'*spdiags(w,0,opt.N,opt.N)*ws;
  data.mass = mass;

function[tf] = mass_validator(data,opt)
  
  mass = data.mass;
  tol = 10^(-6 + opt.s/5);
  if mod(opt.N,2)==0 % Degenerate case
    tf = norm(mass(2:end,2:end)-eye(opt.N-1))<tol;
  else
    tf = norm(mass-eye(opt.N))<tol;
  end

function[data] = interpolation_data(opt)

  import_package('speclab');
  wiener = speclab.wiener;
  %wiener = packages.speclab.wiener;
  [x,w] = wiener.quad.pi_gauss_quadrature(opt.N,opt);
  ks = speclab.common.integer_range(opt.N);
  ws = wiener.eval.wiener_function(x,ks,opt);

  sss = speclab.common.standard_scaleshift_1d;
  f = @(x) exp(-sss(x,opt).^2)./(1+sss(x,opt).^2);
  df = @(x) (-2*sss(x,opt).*f(x) - ...
     2*sss(x,opt).*exp(-sss(x,opt).^2)./(1+sss(x,opt).^2).^2)/opt.scale;

  fx = f(x);
  modes = ws'*(fx.*w);
  
  x_refined = linspace(-3*opt.scale+opt.shift, 3*opt.scale+opt.shift,300).';
  ws_refined = wiener.eval.wiener_function(x_refined,ks,opt);
  fx = f(x_refined);
  [data.modes, data.x_refined, data.ws_refined, data.fx] = deal(modes, ...
     x_refined, ws_refined, fx);

function[tf] = interpolation_validator(data,opt)
  
  tol = 10^(-6+opt.s/5 + abs(opt.scale-1));
  [modes, x_refined, ws_refined, fx] = deal(data.modes,...
  data.x_refined, data.ws_refined, data.fx);

  fx_interp = ws_refined*modes;

  tf = norm(fx_interp-fx)<tol;

function[data] = derivative_data(opt)

  import_package('speclab');
  wiener = speclab.wiener;
  %wiener = packages.speclab.wiener;
  [x,w] = wiener.quad.pi_gauss_quadrature(opt.N,opt);
  ks = speclab.common.integer_range(opt.N);
  ws = wiener.eval.wiener_function(x,ks,opt);

  sss = speclab.common.standard_scaleshift_1d;
  f = @(x) exp(-sss(x,opt).^2)./(1+sss(x,opt).^2);
  df = @(x) (-2*sss(x,opt).*f(x) - ...
     2*sss(x,opt).*exp(-sss(x,opt).^2)./(1+sss(x,opt).^2).^2)/opt.scale;

  fx = f(x);
  modes = ws'*(fx.*w);
  
  x_refined = linspace(-3*opt.scale+opt.shift, 3*opt.scale+opt.shift,300).';
  dws_refined = wiener.eval.derivative_wiener_function(x_refined,ks,opt);
  dfx = df(x_refined);
  [data.modes, data.x_refined, data.dws_refined, data.dfx] = deal(modes, ...
     x_refined, dws_refined, dfx);

function[tf] = derivative_validator(data,opt)
  
  tol = 10^(-4+opt.s/5 + abs(opt.scale-1)/2);
  [modes, x_refined, dws_refined, dfx] = deal(data.modes,...
  data.x_refined, data.dws_refined, data.dfx);

  dfx_interp = dws_refined*modes;

  tf = norm(dfx_interp-dfx)<tol;
