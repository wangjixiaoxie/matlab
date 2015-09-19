function handle=lt_plot_area(X,Y,color, transparency)
%% DEFAULTS
if ~exist('transparency', 'var');
    transparency=0.4;
end


%% 5/12/15 - LT, plotting like matlab's area(), but with defaults

    handle=area(X,Y);
    set(handle,'FaceColor',color,'EdgeColor',color,'LineWidth',2);
    
    % set transparent
    harea_children=get(handle,'Children');
    set(harea_children,'FaceAlpha',transparency);
    
