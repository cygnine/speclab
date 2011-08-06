function[x,w] = gauss_lobatto_quadrature(self,N,varargin)
% [x,w] = gauss_radau_quadrature(self,N, {r1=self.domain.slices{1}.interval(1), 
%                                         r2=self.domain.slices{1}.interval(2)})
%
%     Computes the n-point Gauss-Lobatto quadrature rule associated with the
%     system of orthogonal polynomials. If no optional inputs 'r1' and 'r2'
%     specifying the location of the fixed points are given, they are assumed to
%     be the endpoints of the domain.

persistent inparse gq opoly_evaluate
if isempty(inparse)
  inparse = inputParser();
  inparse.KeepUmatched = false;

  inparse.addParamValue('r1', self.standard_domain.interval(1));
  inparse.addParamValue('r2', self.standard_domain.interval(2));

  from speclab.d1_utils import gauss_quadrature as gq
  from speclab.d1_utils import opoly_evaluate
end
inparse.parse(varargin{:});
opt = inparse.Results;

[a,b] = self.recurrence(0:(N-1));
a = a(1:N);
b = b(1:N);

Q = self;
Q.domain = Q.standard_domain;
Q.normalization = MonicNormalization.instance();
temp = Q.evaluate([opt.r1; opt.r2],[N-1,N-2]);

%temp = self.evaluate([opt.r1; opt.r2],[N-1,N-2],'normalization',MonicNormalization.instance());
% Can't do the above since self.evaluate uses scaling
%p = ones(2);
%temp = p.*opoly_evaluate([opt.r1; opt.r2], a,b, [N-1,N-2], 0);
%% Manually scale to proper normalization, *sigh*
%btemp = sqrt(b(:));
%A = [1; repmat(1, [max([N-1, 0]) 1])];
%btemp = cumprod(btemp.*A);
%temp = temp*diag(btemp([N N-1]));
%% Now the array temp is monic evaluations of p_{N-1} and p_{N-2} at the
%% endpoints of the standard interval.

modif = inv(temp)*[opt.r1*temp(1,1); opt.r2*temp(2,1)];

% Lobatto modification for Jacobi matrix: 1999_gautschi
a(N) = modif(1);
b(N) = modif(2);

[x,w] = gq(a,b);

% Convert back to physical inteval
[x,w] = self.scale_quadrature(x,w);
