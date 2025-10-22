function stats_disp_f(tbl)

% Input to the function
% [p,tbl,stats] = anovan...
% [sta,var]=mes2way

for lineNo = 2:size(tbl,1)-2
    fprintf('\n%s: F(%d,%d) = %0.3f, P = %0.3f \n',...
        tbl{lineNo,1}, tbl{lineNo,3}, tbl{lineNo,11}, tbl{lineNo,6}, tbl{lineNo,7}); 

end