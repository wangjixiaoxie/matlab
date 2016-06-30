function lt_plot_equalyaxis(hfigs, xaxis, yaxis)
%% LT 10/10/15 - input vector of axis handles (e.g. h=figure()) . will take lowest and highest y and apply to all. will also linkaxes.
% xaxis=1;
% yaxis=1;

% ---- equalize y axis
if yaxis==1;
    axis_handles_all=[];
    ylim_all=[];
    for i=1:length(hfigs);
        try
            axis_handles_all=[axis_handles_all get(hfigs(i), 'CurrentAxes')];
        catch err
            axis_handles_all=[axis_handles_all hfigs(i)];
        end
        
        ylim_all=[ylim_all; get(axis_handles_all(end),'Ylim')];
    end
    
    
    ylim_min=min(ylim_all(:,1));
    ylim_max=max(ylim_all(:,2));
    
    for i=1:length(axis_handles_all);
        ax=axis_handles_all(i);
        
        set(ax, 'Ylim', [ylim_min ylim_max]);
    end
    
end

% ----- equalize x axes
if xaxis==1;
    axis_handles_all=[];
    xlim_all=[];
    for i=1:length(hfigs);
        try
            axis_handles_all=[axis_handles_all get(hfigs(i), 'CurrentAxes')];
        catch err
            axis_handles_all=[axis_handles_all hfigs(i)];
        end
        
        xlim_all=[xlim_all; get(axis_handles_all(end),'Xlim')];
    end
    
    
    xlim_min=min(xlim_all(:,1));
    xlim_max=max(xlim_all(:,2));
    
    for i=1:length(axis_handles_all);
        ax=axis_handles_all(i);
        
        set(ax, 'Xlim', [xlim_min xlim_max]);
    end
    
end

tmp=logical([xaxis yaxis]);
string_tmp='xy';
string_tmp=string_tmp(tmp);
linkaxes(axis_handles_all, string_tmp);

