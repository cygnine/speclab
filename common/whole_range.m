function[ks] = whole_range(N)
% [ks] = whole_range(N)
%
%     Wrapper for Matlab's : operator. Takes N ---> [0, 1, 2, ..., N-1]

ks = [0:N].';
