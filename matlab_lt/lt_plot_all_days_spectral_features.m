%% LT 11/19/13
% from analysis of lt_db_seq_func_day_save_rd3gr35, go to the all_days...
% folder and run this to plot summary figures.

function lt_plot_all_days_spectral_features(metrics,curr_dir, phrase, syllables)
% metric= a cell holding strings describing the metrics you want to
% analyze: e.g. {"entropy", "amplitude"}, must correspond to the name used
% in the "all_days..." folder.
% phrase= string indexing your experiment (e.g. randomWN)
% curr_dir = bird folder.

display('plotting over all days')
i_max=length(metrics);

for i=1:i_max; % for each metric
    for j=1:length(syllables);
                
        cd(['/' curr_dir '/all_days_' metrics{i} '_' phrase])
        
        [batch]=lt_write_all_folder_contents_to_batch(syllables(j)); % put all the day .mat files into a batch file
        
        % to plot above data
        %     batch=input('what is name of batch? ','s');
        %     syllable=input('what is the name of the syllable? ','s');
        
        if strcmp(metrics{i},'fraction')==1;
            lt_db_plot_over_experiment(batch,'r', 1, 1,1,1) % with removal of tukey outliers (the 1 at the end is only for 'fraction')
        else
            lt_db_plot_over_experiment(batch,'r', 1, 1,1,0)
        end
        title([metrics{i} ' vs. day; tukey outliers removed. Syllable: ' num2str(syllables(j))])
        % figure; lt_db_plot_over_experiment(batch,'g', 1, 0,1,0) % without removal of tukey outliers
        % title('entropy vs. day; tukey outliers not removed')
                
        saveas(figure(gcf),[metrics{i} '_over_days_of_' num2str(syllables(j))],'fig')
    end
end
