function []=plottrig(ttimes,handle,height)
    axis(handle);
    
    hold on;
    ntrig=size(ttimes);
    htvec=height*ones(ntrig,1);
    plot(ttimes, htvec, 'k^', 'Linewidth', 2);
    