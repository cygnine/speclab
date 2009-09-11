% MATLAB File : NewtonToMonomialDebug.m
%
% * Creation Date : 2009-06-05
%
% * Last Modified : Fri 05 Jun 2009 03:37:53 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Debugging script for NewtonToMonomial.m

%%%%%%%%%%%%%%%%% Test simple polynomial
cd ..
nc = [1;2;3];
nc = [nc,[2;0;-1]];
x = [0;1;2];
x = [x,[-1;1;0]];

% Exact solution:
mc = [1;-1;3];
mc = [mc,[3;0;-1]];

mcCalc = NewtonToMonomial(nc,x);

fprintf('Error is %3.3f\n', norm(mc - mcCalc));

cd debug
