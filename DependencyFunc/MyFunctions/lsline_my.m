function plt = lsline_my(x,y)

[beta,~,~,~,stats] = regress(y,[x*0+1, x]);
% [beta,stats] = robustfit(x, y);

xx = [min(x); max(x)];

hold on; 
plt = plot(xx,beta(1)+beta(2)*xx,'r');

end