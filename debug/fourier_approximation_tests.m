function[container] = fourier_approximation_tests(container,opt)
% [container] = fourier_approximation_tests(container,opt)
%
%     Appends function approximation tests for Fourier Series to the
%     TestContainer container. The Fourier family to use is specified by the
%     struct of parameters opt.

import debug.*

test = ValidationTest('description', 'Mass matrix',...
                      'parameters', opt,...
                      'validator', @mass_validator,...
                      'data_generator', @mass_data);
container = container.append(test);

test = ValidationTest('description', 'Fourier interpolation',...
                      'parameters', opt,...
                      'validator', @interpolation_validator,...
                      'data_generator', @interpolation_data);
container = container.append(test);

test = ValidationTest('description', 'Fourier derivative',...
                      'parameters', opt,...
                      'validator', @derivative_validator,...
                      'data_generator', @derivative_data);
container = container.append(test);

end

function[data] = mass_data(opt)

  global handles;
  fourier = handles.speclab.fourier;
  ks = handles.speclab.common.integer_range(opt.N);
  [x,w] = fourier.quad.gauss_quadrature(opt.N,opt);

  v = fourier.eval.fseries(x,ks,opt);
  mass = v'*spdiags(w,0,opt.N,opt.N)*v;
  [data.v,data.mass] = deal(v,mass);
end

function[tf] = mass_validator(data,opt)
  
  [v,mass] = deal(data.v,data.mass);
  tol = 1e-6;
  if mod(opt.N,2)==0  % degenerate case
    tf = norm(mass(2:end,2:end)-eye(opt.N-1))<tol;
  else
    tf = norm(mass-eye(opt.N))<tol;
  end
end

function[data] = interpolation_data(opt)
  global handles;
  fourier = handles.speclab.fourier;
  ks = handles.speclab.common.integer_range(opt.N);
  sss = handles.speclab.common.standard_scaleshift_1d;

  f = @(x) exp(sin(sss(x,opt)));
  df = @(x) 1/opt.scale*cos(sss(x,opt)).*exp(sin(sss(x,opt)));
  fint = fourier.interval(opt);

  [x,w] = fourier.quad.gauss_quadrature(opt.N,opt);
  x_refined = linspace(fint(1),fint(2),opt.N*3).';
  v = fourier.eval.fseries(x,ks,opt);
  v_refined = fourier.eval.fseries(x_refined,ks,opt);
  modes = v'*spdiags(w,0,opt.N,opt.N)*f(x);
  fx_refined = f(x_refined);

  [data.modes,data.v_refined,data.fx_refined] = deal(modes,v_refined,fx_refined);
end

function[tf] = interpolation_validator(data,opt)

  [modes,v_refined,fx_refined] = deal(data.modes, data.v_refined, data.fx_refined);
  fx_approx = v_refined*modes;

  tol = 10^(-8 + (opt.gamma+opt.delta)/4);
  tf = norm(fx_approx-fx_refined)<tol;
end

function[data] = derivative_data(opt)
  global handles;
  fourier = handles.speclab.fourier;
  ks = handles.speclab.common.integer_range(opt.N);
  sss = handles.speclab.common.standard_scaleshift_1d;

  f = @(x) exp(sin(sss(x,opt)));
  df = @(x) 1/opt.scale*cos(sss(x,opt)).*exp(sin(sss(x,opt)));
  fint = fourier.interval(opt);

  [x,w] = fourier.quad.gauss_quadrature(opt.N,opt);
  x_refined = linspace(fint(1),fint(2),opt.N*3).';
  v = fourier.eval.fseries(x,ks,opt);
  dv_refined = fourier.eval.dfseries(x_refined,ks,opt);
  modes = v'*spdiags(w,0,opt.N,opt.N)*f(x);
  dfx_refined = df(x_refined);

  [data.modes,data.dv_refined,data.dfx_refined] = deal(modes,dv_refined,dfx_refined);
end

function[tf] = derivative_validator(data,opt)

  [modes,dv_refined,dfx_refined] = deal(data.modes, data.dv_refined, data.dfx_refined);
  dfx_approx = dv_refined*modes;

  tol = 10^(-4 + (opt.gamma+opt.delta)/2);
  tf = norm(dfx_approx-dfx_refined)<tol;
end
