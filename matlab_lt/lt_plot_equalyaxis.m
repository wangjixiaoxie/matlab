function lt_plot_equalyaxis(hfigs)
%% LT 10/10/15 - input vector of axis handles (e.g. h=figure()) . will take lowest and highest y and apply to all. will also linkaxes.


% ---- equalize y axis
axis_handles_all=[];
ylim_all=[];
for i=1:length(hfigs);
    try
    axis_handles_all=[axis_handles_all get(hfigs(i), 'CurrentAxes')];
    catch err
    axis_handles_all=[axis_handles_all hfigs(i)];
    end
    
    ylim_all=[ylim_all get(axis_handles_all(end),'Ylim')];
end


    ylim_min=min(ylim_all(:,1));
    ylim_max=max(ylim_all(:,2));
    
    for i=1:length(axis_handles_all);
        ax=axis_handles_all(i);
        
        set(ax, 'Ylim', [ylim_min ylim_max]);
    end
    
    linkaxes(axis_handles_all, 'y');

