function     lt_neural_PLOT_rasterline(spktimes, yval, plotcol, plotdot)

if ~exist('plotdot', 'var')
    plotdot = [];
end

%% plots a single line of raster (lt, 11/4/17)
% spktimes; array, in sec
% yval: position to plot (will go from -1 to ...)



if plotdot==1
   plot(spktimes, yval, '.', 'Color', plotcol); 
else
    for ttt =1:length(spktimes)
        
        line([spktimes(ttt) spktimes(ttt)], [yval-0.4 yval+0.4], ...
            'Color', plotcol, 'LineWidth', 1);
    end
end
