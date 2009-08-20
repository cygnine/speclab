function[alphas,betas] = chebyshev_modified_moments(nus,as,bs);

% [alphas,betas] = chebyshev_modified_moments(nus,as,bs);
% Uses the `Chebyshev Modified Moment' algorithm to derive the recurrence
% quantities for orthogonal polynomials with monomial moments nus.
% The as and bs are the recurrence quantities for some known set of orthogonal
% polynomials under a different weight. 
% The output alphas and betas have length N, where N is the number of modified
% moments nus. 
%
% 20080521: acn

% The number of recurrence quantites we can compute:
N = floor(length(nus)/2);

if (length(as)<2*N) | (length(bs)<2*N);
  error('Need more known recurrence quantities as and/or bs');
end

% Other preprocessing:
as = as(:).';
bs = bs(:).';

alphas = zeros([N 1]);
betas = zeros([N 1]);

sigmas = zeros([2 2*N]);

% Initialization
sigmas(2,:) = 0;
sigmas(1,:) = nus(1:2*N);
sigtemp = zeros([1 2*N]);

alphas(1) = as(1) + nus(2)/nus(1);
betas(1) = nus(1);

for k = 1:(N-1);
  inds = (k+1):(2*N-k);
  sigtemp(:) = 0;
  sigtemp(inds) = sigmas(1,inds+1) - (alphas(k) - as(inds)).*sigmas(1,inds) - ...
  betas(k)*sigmas(2,inds) + bs(inds).*sigmas(1,inds-1);

  sigmas(2,:) = sigmas(1,:);
  sigmas(1,:) = sigtemp(1,:);

  alphas(k+1) = as(k+1) - sigmas(2,k+1)/sigmas(2,k) + ...
  sigmas(1,k+2)/sigmas(1,k+1);

  betas(k+1) = sigmas(1,k+1)/sigmas(2,k);
end

betas(1) = nus(1);
