%% LT 7/6/14 - Modified from 11/18 in 2 ways:
% 1) is now a function
% 2) more drastic changes, so this is v2.  i.e. there is no script version
% of this that is up to day.  e.g. changing how parameters are inputed.
% goal is for this to be incorporated into lt_all_days_various...
% 3) allow syllable to be motif (e.g. ab, not just b) using regular
% expressions
% 4) also removed fvalsstr_all calculation.



%% created by LT 11/18/13 to take individual days and calculate pitch contour and time average pitch
% modified from db_contour...v3, but don't have to do every day.  

% go to the day folder, and do lt_db_transfer_calls(2). Then run this.

% clear all
% close all

function day_pitch=lt_calc_day_pitch_v2_FUNCTION(syl_target, syl_pre, phrase, freq_range, pc_time_window, plot_result,pc_window,varargin);
% EXAMPLES: 
% 1) syl_target = 'ab' (perform for both a and b)
% 2) syl_pre = {'c','d'} (preceding syllables. leave as '' if not needed)
% 2) phrase = 'CPseq' (indexing both data folders and the future save folder)
% 3) freq_range = {[2000 2500],[1000 1500]} - one for each syllable - Important: keep this as tight as possible.
% 4) pc_time_window = {[200 400], [100 300]} i.e. 20 to 40 ms to average over.  If you knwo
% this enter it.  Otherwise enter {0,0} and this function will plot PCs then
% prompt you to enter time windows.
% 5) plot_result = 1 (plot) or 0 (skip)
% 6) pc_window = window for pitch contour - default is 3000.  larger is
% slower, but can see for of the end of the syllable.
% 7) varargin{1}=batch; i.e. types batch name.  if no entry, then will use batch.catch.keep as default.

%% PARAMETERS I WANT TO KEEP PERMANENT

pc_harms=[1]; % harms to weigh in jc_pitch... DB used [1 2] but not worth it I think.
calcFFT=0; % 1 will get fft from findwnote.  keep 0 since that FFT is just at a single point (relteive to onset), waste of time.

if ~isempty(varargin)>0;
    batch=varargin{1};
else
    batch='batch.catch.keep';
end

%% FUNCTION SPECIFIC - converting arguments into parameters originally used int he string version (v1)

if ~exist('pc_window');
pc_window=4000;
end


% [~, ~, date, ~]=lt_get_birdname_date_from_dir(1);

day='08Jul2014';

for i = 1:length(syl_target);
    if strcmp(syl_pre{i},'');
        syllables{i}=syl_target(i);
    elseif strcmp(syl_pre{i},'-');
        syllables{i} = syl_target(i); % dash cannot go into field name.
    else
        syllables{i} = ['pre_' syl_pre{i} '_targ_' syl_target(i)]; % this is the name that combines both pre and target.
    end
    
    frequency_range.(syllables{i}) = freq_range{i};
%     refractory_period.(syllables{i}) = 0.2;  % Dont think this is important
end

%% making save directories (if not existent)
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

%% Calculating 1) fft using findwnote2tw (over a window specified by refractory period, and NFFT), then 2) that fft data is entered into jc_pitchcontourFV to calculate pitch contours.

for i = 1:length(syl_target);
%     fvalsstr_all.(syllables{i}) =findwnote2tw('all_cbin_not_mat',...
%         syllables{i},'',...
%         refractory_period.(syllables{i}),...
%         frequency_range.(syllables{i}),1024,0,'obs0');
    
%     fvalsstr_forpc_all.(syllables{i}) =...
%         findwnote2tw('all_cbin_not_mat',...
%         syllables{i},'',-0.01,...
%         frequency_range.(syllables{i}),8000,0,'obs0');

% try % some days might not have this syllable
    fvalsstr_forpc_all.(syllables{i}) =...
        findwnote2tw_v4_LT(batch,...
        syl_target(i),syl_pre{i},0,...
        frequency_range.(syllables{i}),pc_window,1,1,'obs0',calcFFT); % used to be 8000 but now 3000 as it is much faster in pc step (i.e. how much data to look at)
    
    
    [pc_all.(syllables{i}), pc_F, pc_T] = ...
        jc_pitchcontourFV_LT(fvalsstr_forpc_all.(syllables{i}),...
        1024,1020,1, min(frequency_range.(syllables{i})),...
        max(frequency_range.(syllables{i})),pc_harms,'obs0');
% catch err
% end
end

fprintf('done!\n')


%% Figure of pitch contour comparing mean (and sd) over time

if plot_result==1;
fprintf('\nCreating figures of pitch contours over time...')
for i = 1:length(syllables);
    hfig(i)=figure; hold on
    %     title([nameofbird '   duration: '...
    %         num2str(multiple_pitch.parameters.days{4}-multiple_pitch.parameters.days{3}+1) ' days'...
    %         '   syllable: ' multiple_pitch.parameters.syllables{i}])
    %plots the mean
    plot(pc_all.(syllables{i}),'LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
    plot(mean(pc_all.(syllables{i})'),'Linewidth',2)
   
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
    
    xlabel(['Time, bins: ' num2str(pc_T(2)-pc_T(1))])
    ylabel('Frequency (Hz)')
    title(['Mean and individual pitch contours for ' syllables{i} ', on ' day]); 
    saveas(figure(i), ['../individual_day_pitch_' phrase '/figures/PC_' day '_' syllables{i}], 'fig')
    % plot with real time axis
    hfig_ms(i)=figure; hold on
    plot(pc_T*1000,pc_all.(syllables{i}),'LineStyle','--','Color',[0.6 0.6 0.6]) % plot all pitch contours in light shade
    plot(pc_T*1000,mean(pc_all.(syllables{i})'),'Linewidth',2)
    xlabel('Time (ms)')
    ylabel('Frequency (Hz)')
    title(['Mean and individual pitch contours for ' syllables{i} ', on ' day]); 
    saveas(figure(i), ['../individual_day_pitch_' phrase '/figures/PC_realtime_' day '_' syllables{i}], 'fig')
end

fprintf('done!\n')
end


%% Inputting time range to average pitch contour.  (based on PC figure)
for i = 1:length(syllables)
    if pc_time_window{i}==0;
        time_range.(syllables{i})=input(['time range to average PC over for syllable ' num2str(syllables{i}) ' ?']);
    else
        time_range.(syllables{i})=pc_time_window{i};
        disp('time window/s: '); pc_time_window % show the time window
    end
end

%% Overlaying those time ranges onto PC figures;
if plot_result==1;
for i=1:length(hfig);
    figure(hfig(i));
    line([time_range.(syllables{i})(1) time_range.(syllables{i})(1)], ylim,'Color',[0.5 0.2 0.2]);
    line([time_range.(syllables{i})(2) time_range.(syllables{i})(2)], ylim,'Color',[0.5 0.2 0.2]);
end
end



%% Calculating summary statistics.
% Calculates FF value, time during day, and gets rid of outliers
fprintf('\nCalculating FF over time and getting rid of outliers...')
for i = 1:length(syllables)
    %calculates time of note during the day
    FF.(syllables{i}).time_and_FF(:,1) = ...
        db_timing4(fvalsstr_forpc_all.(syllables{i}))-floor(db_timing4(fvalsstr_forpc_all.(syllables{i})(1)));
        
    %calculates FF of each note during day
    FF.(syllables{i}).time_and_FF(:,2) = ...
        mean(pc_all.(syllables{i})(min(time_range.(syllables{i})):...
        max(time_range.(syllables{i})),:))';
    
    %gets rid of outliers (Tukey's method)
    [FF.(syllables{i}).time_and_FF_NoOutliers, FF.(syllables{i}).outliers_high, FF.(syllables{i}).outliers_low] = ...
        db_tukey_outlier(FF.(syllables{i}).time_and_FF,2);
    
    %calculates mean for day
    FF.(syllables{i}).mean_FF = ...
        mean(FF.(syllables{i}).time_and_FF(:,2));
    
    % median
    FF.(syllables{i}).median_FF = ...
        median(FF.(syllables{i}).time_and_FF(:,2));
   
    %calculates sd for day
    FF.(syllables{i}).sd_FF = ...
        std(FF.(syllables{i}).time_and_FF(:,2));
    
    % calculates 70th percentiles:
    FF.(syllables{i}).prctiles_5_30_70_95=prctile(FF.(syllables{i}).time_and_FF(:,2),[5 30 70 95]);
    
end
fprintf('done!\n')

%% Scatter plot of pitch over day
if plot_result==1;
fprintf('\nCreating figures of FF over time...')

for i = 1:length(syllables)
    figure(), hold on
    title([day '; syllable: ' syllables{i} '. Time window: ' num2str(time_range.(syllables{i}))])
    ylabel('Frequency (Hz)')
    xlabel('syllable rendition #')
    
    %     if isempty(multiple_pitch.FF.(multiple_pitch.parameters.syllables{i}).time_and_FF{end}) == 0
    %         xlim([-1-multiple_pitch.parameters.duration.begin...
    %             ceil(multiple_pitch.FF.(multiple_pitch.parameters.syllables{i}).time_and_FF{end}(end,1))+1])
    %     else
    %         xlim([-1-multiple_pitch.parameters.duration.begin ...
    %            -1-multiple_pitch.parameters.duration.begin + length(multiple_pitch.FF.(multiple_pitch.parameters.syllables{i}).time_and_FF)])
    %     end
    
    %     %plots FF for time range
    %     plot(FF.(syllables{i}).time_and_FF(:,1),...
    %         FF.(syllables{i}).time_and_FF(:,2),...
    %         '.')

    % plots with x-values being 1,2,3..., instead of rendition time.
    plot(FF.(syllables{i}).time_and_FF(:,2),...
        '.')
    
    %plots FF mean per day
    plot(median(FF.(syllables{i}).time_and_FF(:,1)),...
        FF.(syllables{i}).mean_FF,...
        'ko',...
        'MarkerFaceColor','k',...
        'LineWidth',2)
    
    %plots FF sd per day
    errorbar(median(FF.(syllables{i}).time_and_FF(:,1)),...
        FF.(syllables{i}).mean_FF,...
        FF.(syllables{i}).sd_FF,...
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
    %             [multiple_pitch.FF.(multiple_pitch.parameters.syllables{i}).baseline_mean...
    %             multiple_pitch.FF.(multiple_pitch.parameters.syllables{i}).baseline_mean],...
    %             'Color', [.3 .3 .3],...
    %             'LineStyle', '--')
    %     catch err
    %         continue
    %     end
    
    saveas(figure(gcf), ['../individual_day_pitch_' phrase '/figures/FF_over_time_' day '_' syllables{i}], 'fig')
    cd(current_directory)
end
fprintf('done!\n')
end
%% assigning structure fields to all important variables

day_pitch.parameters.syllables=syllables;
day_pitch.parameters.day=day;
day_pitch.parameters.phrase=phrase;
day_pitch.parameters.frequency_range=frequency_range;
% day_pitch.parameters.refractory_period=refractory_period;
day_pitch.parameters.time_range=time_range;
day_pitch.parameters.current_directory=current_directory;
day_pitch.parameters.syl_pre=syl_pre;
day_pitch.parameters.syl_target=syl_target;


% day_pitch.fvalsstr_all=fvalsstr_all;
day_pitch.fvalsstr_forpc_all=fvalsstr_forpc_all;
day_pitch.pc_all=pc_all;
day_pitch.FF=FF;


%% Saving
cd(['../individual_day_pitch_' phrase])

save(['day_pitch_' day],'day_pitch')


cd(current_directory)

%% Printing stats:
day
disp('mean FF');
disp('STD of FF');
disp('sample size');
for i = 1:length(syllables);
    disp(['syllable: ' num2str(syllables{i})])
    disp(day_pitch.FF.(syllables{i}).mean_FF)
    disp(day_pitch.FF.(syllables{i}).sd_FF)
    disp(100*(day_pitch.FF.(syllables{i}).sd_FF)/(day_pitch.FF.(syllables{i}).mean_FF))
    disp(length(day_pitch.FF.(syllables{i}).time_and_FF))
    disp(FF.(syllables{i}).prctiles_5_30_70_95)
end
    