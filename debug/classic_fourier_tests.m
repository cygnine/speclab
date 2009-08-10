function[container] = classic_fourier_tests(container,opt);
% [container] = classic_fourier_tests(container,opt);
%
%     Appends some Fourier Series tests for the classical Fourier Series. The
%     Fourier family is defined by opt, and container is of type TestContainer.

import debug.*

test = ValidationTest('description', 'Evaluation same as classic FS',...
                      'parameters', opt,...
                      'validator', @eval_fs_validator,...
                      'data_generator', @eval_fs_data);
container = container.append(test);

test = ValidationTest('description', 'Derivative same as classic FS',...
                      'parameters', opt,...
                      'validator', @derivative_fs_validator,...
                      'data_generator', @derivative_fs_data);
container = container.append(test);

end

function[data] = eval_fs_data(opt)
  global handles;
  fourier = handles.speclab.fourier;
  sss = handles.speclab.common.standard_scaleshift_1d;

  fint = fourier.interval(opt);
  x = linspace(fint(1),fint(2),300).';
  ks = handles.speclab.common.integer_range(opt.N);
  ks = reshape(ks,[1,length(ks)]);

  fs = fourier.eval.fseries(x,ks,opt);
  fs_exp = 1/sqrt(2*pi)*exp(i*sss(x,opt)*ks);
  [data.fs, data.fs_exp] = deal(fs,fs_exp);
end

function[tf] = eval_fs_validator(data,opt)
  tol = 1e-8;
  [fs,fs_exp] = deal(data.fs, data.fs_exp);

  tf = norm(fs - fs_exp)<tol;
end

function[data] = derivative_fs_data(opt)
  global handles;
  fourier = handles.speclab.fourier;
  sss = handles.speclab.common.standard_scaleshift_1d;

  fint = fourier.interval(opt);
  x = linspace(fint(1),fint(2),300).';
  ks = handles.speclab.common.integer_range(opt.N);
  ks = reshape(ks,[1,length(ks)]);

  dfs = fourier.eval.dfseries(x,ks,opt);
  dfs_exp = 1/sqrt(2*pi)*exp(i*sss(x,opt)*ks)*spdiags(i*ks.',0,opt.N,opt.N);
  dfs_exp = dfs_exp/opt.scale;
  [data.dfs, data.dfs_exp] = deal(dfs,dfs_exp);
end

function[tf] = derivative_fs_validator(data,opt)
  tol = 1e-8;
  [dfs,dfs_exp] = deal(data.dfs, data.dfs_exp);

  tf = norm(dfs - dfs_exp)<tol;
end
