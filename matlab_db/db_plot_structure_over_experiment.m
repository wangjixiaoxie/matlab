function [] = db_plot_structure_over_experiment( batch, syl, color, hold_off, tukey, plot_all )
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

if isempty(findobj('type','figure')) == 1
    figure()
else
end

hold on

hold on

fid = fopen(batch,'r');

tline = fgetl(fid);
i = 1;

while ischar(tline)
    load(tline)
    
    if i == 1
        if strcmpi(syl,'');
            syl = fieldnames(data);
        else
            syl = {syl};
        end
        
        if isempty(time.(syl{i})) == 1
            find_dashes = strfind(tline,'_');
            start_date = tline(find_dashes(1)+1:find_dashes(2)-1);
            start_date = datenum(start_date)-1;
        else
            start_date = min(floor(time.(syl{i})))-1;
        end
    end
    

    for j = 1:max(size(syl))
        if isempty(data.(syl{j})) == 1
        else
            if tukey == 1;
                %gets rid of tukey outliers (3.5*IQR)
                [tukey_data, high, low] = db_tukey_outlier(data.(syl{j}));
                tukey_time = removerows(time.(syl{j}), 'ind', [high;low]);
                
                if plot_all == 1
                    %plots data points
                    plot(tukey_time-start_date,tukey_data,'o','Color',color(j))
                end
                hold on
                
                %plots median of data
                plot(median(tukey_time)-start_date, mean(tukey_data), 'o',...
                    'Color', color(j), 'MarkerFaceColor', color(j), 'MarkerSize', 10, 'MarkerEdgeColor', 'k')
                
                %plots std errorbars
                errorbar(median(tukey_time)-start_date, mean(tukey_data), std(tukey_data),...
                    'Color', color(j), 'LineWidth', 2)
            else
                if plot_all == 1
                    plot(time-start_date,data,'o','Color',color(j))
                end
                
                %plots median of data
                plot(median(time)-start_date, mean(data), 'o',...
                    'Color', color(j), 'MarkerFaceColor', color(j), 'MarkerSize', 10,'MarkerEdgeColor', 'k' )
                
                %plots std errorbars
                errorbar(median(time)-start_date, mean(data), std(data),...
                    'Color', color(j), 'LineWidth', 2)
            end
        end
    end

    tline = fgetl(fid);
    i = i+1;
end

if hold_off == 1
    hold off
end


end

