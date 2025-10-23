function [reverseStr] = display_status(process, reverseStr, current, total)
% Display status
msg = sprintf('%s: %d/%d\n', process, current, total);
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));
end