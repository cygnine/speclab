function[] = disp(self)
% disp -- Display method for Normalization
%
% disp(self)
%
%     Displays description and id information for this instance.

fprintf(['\n' self.description, '\n\nThe following strings are equivalent to this class via the equality (==) operator:\n']);
fprintf('  %s\n', self.ids{:});
fprintf('\n');

