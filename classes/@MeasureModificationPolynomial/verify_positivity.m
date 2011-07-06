function[self] = verify_positivity(self)
% verify_positivity -- Verifies and classifies roots to ensure positivty
%
% self = verify_positivity(self)
%
%     This method should be called whenever self.roots is updated. This method
%     (a) categorizes the roots in self.roots into exterior,
%     (complex-)conjugate, and quadratic, (b) ensures the proper counting of
%     those roots, and (c) updates self.leading_coefficient to ensure positivity
%     on self.interval.

% Make sure roots define a positive polynomial on self.interval, and then
% categorize roots

% Eliminate imag parts less then self.tol
roots = self.roots;
realrootids = abs(imag(roots)) < self.tol;
roots(realrootids) = real(roots(realrootids));

% Now do initial categorization of roots
exteriorids = realrootids | (real(roots) >= self.interval.interval(2)) | ...
                            (real(roots) <= self.interval.interval(1));
ccids = not(realrootids);
quadids = not(exteriorids | ccids);

% By definition, exterior roots are verified
self.exterior_roots = roots(exteriorids);

% Quadratic roots must come in pairs
quadroots = roots(quadids);
if mod(numel(quadroots), 2) == 1
  error('Roots inside in the interval do not come in quadratic pairs -- polynomial is not positive on the interval');
end
self.quadratic_roots = zeros([numel(quadroots)/2 1]);
rootcount = 1;
while numel(quadroots) > 0
  currroot = quadroots(1);
  matchingrootid = find(abs(quadroots(2:end) - currroot) < self.tol, 1, 'first');
  if isempty(matchingrootid)
    error('Roots inside in the interval do not come in quadratic pairs -- polynomial is not positive on the interval');
  end
  self.quadratic_roots(rootcount) = currroot;
  quadroots([1 matchingrootid]) = [];
  rootcount = rootcount + 1;
end
roots(quadids) = reshape(repmat(self.quadratic_roots.', [2 1]), [2*(rootcount-1) 1]);

% complex roots must be conjugate
ccroots = roots(ccids);
if mod(numel(ccroots), 2) == 1
  error('Complex roots do not come in c.c. pairs -- polynomial is not positive on the interval');
end
self.conjugate_roots = zeros([numel(ccroots)/2 1]);
rootcount = 1;
while numel(ccroots) > 0
  currroot = ccroots(1);
  matchingrootid = find(abs(ccroots(2:end) - currroot) < self.tol, 1, 'first');
  if isempty(matchingrootid)
    error('Complex roots do not come in c.c. pairs -- polynomial is not positive on the interval');
  end
  self.conjugate_roots(rootcount) = currroot;
  ccroots([1 matchingrootid]) = [];
  rootcount = rootcount + 1;
end
roots(ccids) = reshape([self.conjugate_roots.' self.conjugate_roots'], [2*(rootcount-1) 1]);

% Update things
self.roots = roots;
self.coefficients = real(self.coefficients);

pval = self.evaluate(self.interval.centroid);
if pval < 0
  warning('Given leading coefficient makes polynomial negative on interval -- changing the sign');
  self.leading_coefficient = -self.leading_coefficient;
end
