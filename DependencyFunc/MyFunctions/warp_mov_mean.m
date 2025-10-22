function k = warp_mov_mean(y,win)
    y_mod = [y(end-1), y(end),y,y(1),y(2)];
  
    k = movmean(y_mod,win);
    k([1,2,end-1,end]) = [];

return