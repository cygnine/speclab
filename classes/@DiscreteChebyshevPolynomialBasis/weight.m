function[w] = weight(self, x)
% weight -- Weight function for Discrete Chebyshev polynomials
%
% w = weight(self, x)
% 
%     On the standard interval [0,N-1], define the weight function
%
%       w(r) = \sum_{n=0}^{N-1} \delta_n(r)
%
%     On the physical domain x(r), this function evaluates C*w(r(x)), where C is
%     determined by the scaling in self.scale_weight and
%     self.weight_normalization. Note that \delta_n is the Dirac delta function
%     collocated at r=n. Since we don't want Inf's floating around, we just
%     return constant of proportionality in front of the Dirac measure at r=n.

r = self.map_to_standard_domain(x);
w = (r==floor(r));
w(r<0 | r > self.N-1) = 0;
% Dirac measures: scaling the input space changes normalization of output space
w = w*self.map_to_domain.A;
w = self.scale_weight(w);
