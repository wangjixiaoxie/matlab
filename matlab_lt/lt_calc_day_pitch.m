%% NOTE: RECOMMENDED to use v2_FUNCTION, since there are many improvements there (not just change to function)

%% created by LT 11/18/13 to take individual days and calculate pitch contour and time average pitch
% modified from db_contour...v3, but don't have to do every day.  

% go to the day folder, and do lt_db_transfer_calls(2)

% clear all
% close all

%% inputting variables
syllables=input(['What are the syllables (e.g. ab)? '], 's');
[birdname, bluejaynum, date, phrase]=lt_get_birdname_date_from_dir(1);
day=date{2};

for i = 1:length(syllables);
    frequency_range.(syllables(i)) =...
        input(['Frequency range for ' syllables(i) '?\n(format [start stop])  ']);
    refractory_period.(syllables(i)) = 0.01;  % input(['Refractory period for ' syllables(i) ' for?  ']);
end

%% output data to:
% fvalsstr_all;
% fvalsstr_forpc_all=[]
% pc_all=[]

%% making save directories
current_directory=pwd;
cd ..
try
    cd(['individual_day_pitch_' phrase])
catch err
    mkdir(['individual_day_pitch_' phrase]);
    cd(['individual_day_pitch_' phrase]);
end


try
    cd('figures')
catch err
    mkdir('figures');
end

cd(current_directory)

%% Calculating pitch contours

% this part takes a long time to process (specifically, the jc_pitchcontour
% analysis)


for i = 1:length(syllables);
    fvalsstr_all.(syllables(i)) =findwnote2tw('all_cbin_not_mat',...
        syllables(i),'',...
        refractory_period.(syllables(i)),...
        frequency_range.(syllables(i)),1024,0,'obs0');
    
    fvalsstr_forpc_all.(syllables(i)) =...
        findwnote2tw('all_cbin_not_mat',...
        syllables(i),'',-0.01,...
        frequency_range.(syllables(i)),8000,0,'obs0');
    
    pc_all.(syllables(i)) = ...
        jc_pitchcontourFV(fvalsstr_forpc_all.(syllables(i)),...
        1024,1020,1, min(frequency_range.(syllables(i))),...
        max(frequency_range.(syllables(i))),[1 2],'obs0');
end

fprintf('done!\n')

%% Figure of pitch contour comparing mean (and sd) over time

fprintf('\nCreating figures of pitch contours over time...')

for i = 1:length(syllables)
    figure(i), hold on
%     title([nameofbird '   duration: '...
%         num2str(multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}+1) ' days'...
%         '   syllable: ' multiple_pitch.parameters.syllables(i)])
    xlabel('Time (10^-4 sec)')
    ylabel('Frequency (Hz)')
        %plots the mean
            plot(mean(pc_all.(syllables(i))'),'Linewidth',2)
            %             %plots mean + sd
            %             plot(mean(pc_all.(nameofbirds{jj}).(syllables{jj}(i)){j-days.(nameofbirds{jj}){3}+1}')...
            %                 +std(pc_all.(nameofbirds{jj}).(syllables{jj}(i)){j-days.(nameofbirds{jj}){3}+1}'),...
            %                 'Color',[1-(j-days.(nameofbirds{jj}){3})./(days.(nameofbirds{jj}){4} - days.(nameofbirds{jj}){3})...
            %                 0 0+(j-days.(nameofbirds{jj}){3})./(days.(nameofbirds{jj}){4} - days.(nameofbirds{jj}){3})])
            %             %plots mean -sd
            %              plot(mean(pc_all.(nameofbirds{jj}).(syllables{jj}(i)){j-days.(nameofbirds{jj}){3}+1}')...
            %                 -std(pc_all.(nameofbirds{jj}).(syllables{jj}(i)){j-days.(nameofbirds{jj}){3}+1}'),...
            %                 'Color',[1-(j-days.(nameofbirds{jj}){3})./(days.(nameofbirds{jj}){4} - days.(nameofbirds{jj}){3})...
            %                 0 0+(j-days.(nameofbirds{jj}){3})./(days.(nameofbirds{jj}){4} - days.(nameofbirds{jj}){3})])
    
        saveas(figure(i), ['../individual_day_pitch_' phrase '/figures/PC_' day '_' syllables(i)], 'fig')
end

fprintf('done!\n')

%% inputting time range to average pitch contour.  (based on PC figure)
for i = 1:length(syllables)
    time_range.(syllables(i))=input(['time range to average PC over for syllable ' num2str(syllables(i)) ' ?']);
end


%% averaging pitch contrours over time.
% Calculates FF value, time during day, and gets rid of outliers
fprintf('\nCalculating FF over time and getting rid of outliers...')
for i = 1:length(syllables)
    %calculates time of note during the day
    FF.(syllables(i)).time_and_FF(:,1) = ...
        db_timing4(fvalsstr_all.(syllables(i)))-floor(db_timing4(fvalsstr_all.(syllables(i))(1)));
        
    %calculates FF of each note during day
    FF.(syllables(i)).time_and_FF(:,2) = ...
        mean(pc_all.(syllables(i))(min(time_range.(syllables(i))):...
        max(time_range.(syllables(i))),:))';
    
%     %gets rid of outliers (Tukey's method)
%     [FF.(syllables(i)).time_and_FF, FF.(syllables(i)).outliers_high, FF.(syllables(i)).outliers_low] = ...
%         db_tukey_outlier(FF.(syllables(i)).time_and_FF,2);
    
    %calculates mean for day
    FF.(syllables(i)).mean_FF = ...
        mean(FF.(syllables(i)).time_and_FF(:,2));
    
    %calculates sd for day
    FF.(syllables(i)).sd_FF = ...
        std(FF.(syllables(i)).time_and_FF(:,2));
    
end
fprintf('done!\n')

%% figure

fprintf('\nCreating figures of FF over time...')

for i = 1:length(syllables)
    figure(), hold on
    title([day '; syllable: ' syllables(i)])
    ylabel('Frequency (Hz)')
    xlabel('syllable rendition #')
    %     if isempty(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{end}) == 0
    %         xlim([-1-multiple_pitch.parameters.duration.begin...
    %             ceil(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF{end}(end,1))+1])
    %     else
    %         xlim([-1-multiple_pitch.parameters.duration.begin ...
    %            -1-multiple_pitch.parameters.duration.begin + length(multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).time_and_FF)])
    %     end
    
    %     %plots FF for time range
    %     plot(FF.(syllables(i)).time_and_FF(:,1),...
    %         FF.(syllables(i)).time_and_FF(:,2),...
    %         '.')
    %
    % plots with x-values being 1,2,3..., instead of rendition time.
    plot(FF.(syllables(i)).time_and_FF(:,2),...
        '.')
    
    %plots FF mean per day
    plot(median(FF.(syllables(i)).time_and_FF(:,1)),...
        FF.(syllables(i)).mean_FF,...
        'ko',...
        'MarkerFaceColor','k',...
        'LineWidth',2)
    
    %plots FF sd per day
    errorbar(median(FF.(syllables(i)).time_and_FF(:,1)),...
        FF.(syllables(i)).mean_FF,...
        FF.(syllables(i)).sd_FF,...
        'ko',...
        'LineWidth',2)
    
    %plots lines for WN begin, consolidation begin, and WN end
    %line for WN begin
    %     try
    %         line([0-.25 0-.25],...
    %             get(gca, 'YLim'),...
    %             'Color', [0 0 0],...
    %             'LineStyle', '--')
    %     catch err
    %         continue
    %     end
    %     %line for consolidation begin
    %     try
    %         line([multiple_pitch.parameters.duration.shift-.25...
    %             multiple_pitch.parameters.duration.shift-.25],...
    %             get(gca, 'YLim'),...
    %             'Color', [0 0 0],...
    %             'LineStyle', '--')
    %     catch err
    %         continue
    %     end
    %     %line for WN end
    %     try
    %         line([multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation-.25...
    %             multiple_pitch.parameters.duration.shift+multiple_pitch.parameters.duration.consolidation-.25],...
    %             get(gca, 'YLim'),...
    %             'Color', [0 0 0],...
    %             'LineStyle', '--')
    %     catch err
    %         continue
    %     end
    %     %line for baseline pitch
    %     try
    %         line(get(gca, 'XLim'),...
    %             [multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).baseline_mean...
    %             multiple_pitch.FF.(multiple_pitch.parameters.syllables(i)).baseline_mean],...
    %             'Color', [.3 .3 .3],...
    %             'LineStyle', '--')
    %     catch err
    %         continue
    %     end
    
    
    
    saveas(figure(gcf), ['../individual_day_pitch_' phrase '/figures/FF_over_time_' day '_' syllables(i)], 'fig')
    cd(current_directory)
end
fprintf('done!\n')

%% assigning structure fields to all important variables

day_pitch.parameters.syllables=syllables;
day_pitch.parameters.day=day;
day_pitch.parameters.phrase=phrase;
day_pitch.parameters.frequency_range=frequency_range;
day_pitch.parameters.refractory_period=refractory_period;
day_pitch.parameters.time_range=time_range;
day_pitch.parameters.current_directory=current_directory;

day_pitch.fvalsstr_all=fvalsstr_all;
day_pitch.fvalsstr_forpc_all=fvalsstr_forpc_all;
day_pitch.pc_all=pc_all;

day_pitch.FF=FF;


%% saving
cd(['../individual_day_pitch_' phrase])

save(['day_pitch_' day],'day_pitch')


cd(current_directory)

%% printing stats:

disp('mean FF');
disp('STD of FF');
disp('COV');
disp('sample size');
for i = 1:length(syllables);
    disp(['syllable: ' num2str(syllables(i))])
    disp(day_pitch.FF.(syllables(i)).mean_FF)
    disp(day_pitch.FF.(syllables(i)).sd_FF)
    disp(100*(day_pitch.FF.(syllables(i)).sd_FF)/(day_pitch.FF.(syllables(i)).mean_FF))
    disp(length(day_pitch.FF.(syllables(i)).time_and_FF))
end
    