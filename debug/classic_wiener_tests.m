function[container] = classic_wiener_tests(container,opt)
% [container] = classic_wiener_tests(container,opt)
%
%     Appends Wiener tests for the classical (s=1) Wiener basis.

import debug.*

test = ValidationTest('description', 'Evaluation same as classic WF',...
                      'parameters', opt,...
                      'validator', @eval_wiener_validator,...
                      'data_generator', @eval_wiener_data);
container = container.append(test);

function[data] = eval_wiener_data(opt)
  global packages;
  wiener = packages.speclab.wiener;
  sss = packages.speclab.common.standard_scaleshift_1d;

  x = linspace(-10,10,300).';
  ks= packages.speclab.common.integer_range(opt.N);
  ks = reshape(ks,[1,length(ks)]);

  wf = wiener.eval.wiener_function(x,ks,opt);
  wf_form = zeros(size(wf));
  x = sss(x,opt);
  for q = 1:length(ks)
    k = ks(q);
    wf_form(:,q) = 1/sqrt(pi*opt.scale)*((i-x)./(i+x)).^k;
    wf_form(:,q) = wf_form(:,q)./(x-i);
  end
  [data.wf, data.wf_form] = deal(wf,wf_form);

function[tf] = eval_wiener_validator(data,opt)
  tol = 1e-8;
  [wf,wf_form] = deal(data.wf, data.wf_form);

  tf = norm(wf - wf_form)<tol;
