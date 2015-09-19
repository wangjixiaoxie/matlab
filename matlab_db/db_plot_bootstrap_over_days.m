function [] = db_plot_bootstrap_over_days( batchfile , color )
%db_plot_fraction_over_days Plot bootstrap data over days (median + 95% CI)
%   Looks for the variable bootstrap (unless boot_norm = 1, then uses boot_norm)
%   in each loaded file then takes the 5,
%   50, and 95th percentiles and plots them for each day. 

if isempty(findobj('type','figure')) == 1
    figure()
else
end

hold on

fid = fopen(batchfile,'r');

tline = fgetl(fid);
i = 1;

while ischar(tline)
    load(tline)
    

    [percentiles] = [prctile(bootstrap,5) prctile(bootstrap,50) prctile(bootstrap,95)];

    
    
    %plots 95% confidence interval
    line([i i],[percentiles(1) percentiles(3)], 'Color', color, 'LineWidth',2)
    
    %puts side lines on edges of CI
    line([i-.2 i+.2],[percentiles(1) percentiles(1)], 'Color', color, 'LineWidth',2)
    line([i-.2 i+.2],[percentiles(3) percentiles(3)], 'Color', color, 'LineWidth',2)
    
    %plots median
    line([i-.1 i+.1],[percentiles(2) percentiles(2)], 'Color', color, 'LineWidth',4)
    
    i = i+1;
    tline = fgetl(fid);
end
fclose(fid);
    


end

