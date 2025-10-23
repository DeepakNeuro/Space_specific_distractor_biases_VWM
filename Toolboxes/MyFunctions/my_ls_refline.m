function my_ls_refline(x,y)

p = polyfit(x,y,1);
f = polyval(p,x);
plot(x,f,'-r', 'LineWidth', 2)

end