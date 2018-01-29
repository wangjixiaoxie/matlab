function lt_save_all_figs
%% LT 4/7/15 - saves all open figures as figure1, figure2,...
% note, order might not necessarily be correct, and will start from 1, ...


h = get(0,'children');
for i=1:length(h)
    
    if (0)
        %     set(h(i),'units','normalized','outerposition',[0 0 1 1]); % make fullscreen before save
        set(h(i), 'Position', [1081 840 1920 1080])
    end
    saveas(h(i), ['figure' num2str(length(h)+1-i)], 'fig');
    %   saveas(h(i), ['figure' num2str(length(h)+1-i)], 'jpg');
    if (0)
        set(h(i),'units','normalized','outerposition',[0 0 0.5 0.5]); % make smaller
    end
end


