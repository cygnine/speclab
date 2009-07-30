function[flags,descriptions,parameters] = concatenate_debug_information(flags,descriptions,parameters,f,d,p);
% [FLAGS,DESCRIPTIONS,PARAMETERS] = CONCATENATE_DEBUG_INFORMATION(FLAGS,DESCRIPTIONS,PARAMETERS,F,D,P);
%
%     Concatenates the boolean array FLAGS with input F.
%     Concatenates the cell array DESCRIPTIONS with input D.
%     Concatenates the cell array PARAMETERS with input P.

M = length(f);
flags(end+1:end+M) = f;

M = length(d);
descriptions{end+1:end+M} = d;

M = length(p);
parameters{end+1:end+M} = p;
