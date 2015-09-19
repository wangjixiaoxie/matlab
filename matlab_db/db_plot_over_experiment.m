function [] = db_plot_over_experiment( batch, color, hold_off, tukey, plot_all )
%db_plot_over_experiment Takes a type of data (ex: motif duration) and
%plots it over time can use Tukey to get rid of outliers
%   Will go through a batch file of a certain type of data. Variables must
%   be called 'data' and 'time'. Will plot the data by day with mean + sd
%   errorbars. Can select color and to take out outliers (set tukey to 1).
%   To hold off, seet hold_off to 1. Plot_all plots all data points (1), if
%   off (0), plots only mean and sd.

if isempty(findobj('type','figure')) == 1
    figure()
else
end

hold on

fid = fopen(batch,'r');

tline = fgetl(fid);
i = 1;

while ischar(tline)
    load(tline)
    
    if i == 1
        if isempty(time) == 1
            find_dashes = strfind(tline,'_');
            start_date = tline(find_dashes(1)+1:find_dashes(2)-1);
            start_date = datenum(start_date)-1;
        else
            start_date = min(floor(time))-1;
        end
    end
    
    if isempty(data) == 1
    else
        if tukey == 1;
            %gets rid of tukey outliers (3.5*IQR)
            [tukey_data, high, low] = db_tukey_outlier(data);
            tukey_time = removerows(time, 'ind', [high;low]);
            
            if plot_all == 1
                %plots data points
                plot(tukey_time-start_date,tukey_data,'o','Color',color)
            end
            
            %plots median of data
            plot(median(tukey_time)-start_date, mean(tukey_data), 'o',...
                'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerSize', 10)
            
            %plots std errorbars
            errorbar(median(tukey_time)-start_date, mean(tukey_data), std(tukey_data),...
                'Color', [0 0 0], 'LineWidth', 2)
        else
            if plot_all == 1
                plot(time-start_date,data,'o','Color',color)
            end
            
            %plots median of data
            plot(median(time)-start_date, mean(data), 'o',...
                'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerSize', 10)
            
            %plots std errorbars
            errorbar(median(time)-start_date, mean(data), std(data),...
                'Color', [0 0 0], 'LineWidth', 2)
        end
    end
    
    tline = fgetl(fid);
    i = i+1;
end

if hold_off == 1
    hold off
end


end

