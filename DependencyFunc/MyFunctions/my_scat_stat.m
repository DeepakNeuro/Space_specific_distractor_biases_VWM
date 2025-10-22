function [x_sorted, y_mean_sorted, y_sem_sorted]  = my_scat_stat(x, y, radius)

% x and y should both be vectors of equal size. 


y_mean = scatstat1(x,y,radius); 

func_sem = @(xx) std(xx)./sqrt(length(xx)); 
y_sem = scatstat1(x,y,radius, func_sem); 

[~, idx] = sort(x); 

x_sorted = x(idx); 
y_mean_sorted = y_mean(idx);
y_sem_sorted = y_sem(idx);

end