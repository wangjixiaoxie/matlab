%% LT 12/9/13 - to plot data compiled previously in lt_all_days_all_analysis, which compiles all analyses into one strucure.
%% LT 3/14/14 - added all_trans plotting

%% LT 5/17/14 - PLOTTING, for WN learning (any type), hits/song, hits/syl sang, and song lengths, across learning
% use all_days_all_analysis structure. Run this in the bird folder.
%   gives plots of 1) how effective was your template over days? 2) how well
%   did the bird learn to avoid WN (per song and per syl sang).
%   NEEDS: lt_get_all_trans... and triglabel data.

lt_all_days_PLOT_WN_stats('19May2014_2256')


%% PLOT entropy over all transitions and over the desired transition

clear all
load all_days_all_analysis_02Apr2014_1416_NIFlesion_AllTrans_WNtrigLabel

% Input parameters
baseline_days=[1 2 3]; %in 1,2,3,... indices, not days with gaps.
WN_days= 9;
thresh_prob=[0.005 0.01 0.02]; % ignore transitions with prob lower than this number
    
% Figure out the syllables that were sang every day.
for j=1:length(thresh_prob);
    for i=1:size(all_days_all_analysis.data,2);
        try
            dummy=all_days_all_analysis.data{i}.all_songs.all_trans.summary.matrix_of_fraction_of_all_trans(all_days_all_analysis.data{i}.all_songs.all_trans.summary.matrix_of_fraction_of_all_trans>thresh_prob(j));
            entropy.all_FractionOfAll(i,j)=sum(-dummy.*log2(dummy));
            %         below, do for specific transitions (div or conv specific)
            %         dummy2=all_days_all_analysis.data{i}.all_songs.all_trans.summary.divergence.matrix_of_
            %         entropy.all_Conv
        catch err
            entropy.all_FractionOfAll(i,j)=nan;
            continue
        end
    end
end

figure;hold on;
for j=1:length(thresh_prob)
    plot(entropy.all_FractionOfAll(:,j),'Color',[(j/length(thresh_prob))*1 0.5 0.5],'Marker','o');
    legend
end

all_days_all_analysis
line([baseline_days(end)+0.5 baseline_days(end)+0.5],ylim,'Color','g')
line([baseline_days(end)+WN_days+0.5 baseline_days(end)+WN_days+0.5],ylim,'Color','r')

    
    


% entropy


%% TRANSITIONS (from lt_get_all_tra... 

% PLOT transition matrices and difference matrices for all days
con_div = 'divergence';
all_days_all_analysis=lt_all_days_PLOT_AllTrans([1 2 3], '25Apr2014_1343', con_div);
lt_all_days_PLOT_AllTrans_specific_trans(con_div,1,'','15May2014_1809') % run this in bird folder.  

con_div = 'convergence';
all_days_all_analysis=lt_all_days_PLOT_AllTrans([1 2 3], '25Apr2014_1343', con_div);
lt_all_days_PLOT_AllTrans_specific_trans(con_div,1,'','15May2014_1809') % run this in bird folder.  



%% THIRD, plot average matrices for pre, during, and post WN.



%% PLOT OVER DAYS, any type of analysis (except pitch)
close all
analysis_type='gap_duration';
first_day_num=datenum(all_days_all_analysis.parameters.first_day);

if strcmp(analysis_type,'gap_duration')==1;
    transitions=all_days_all_analysis.parameters.gap_duration.transitions;
    
    for j=1:length(transitions);
        figure(j); hold on;
        for i=1:length(all_days_all_analysis.data); % number of days
            try
                num_of_renditions=length(all_days_all_analysis.data{i}.all_songs.(transitions{j}).gap_duration);
                scatter(i-1.25+[1:num_of_renditions]./(2*num_of_renditions),all_days_all_analysis.data{i}.all_songs.(transitions{j}).(analysis_type),'b');
            catch err
                continue
            end
        end
        
        % Put in medians, std
        for i=1:length(all_days_all_analysis.data); % number of days
            try
                errorbar(i-1,median(all_days_all_analysis.data{i}.all_songs.(transitions{j}).(analysis_type)),...
                    std(all_days_all_analysis.data{i}.all_songs.(transitions{j}).(analysis_type)),'Marker','.','MarkerSize',27,'Color','k');
            catch err
                continue
            end
        end
        
        % put in baseline lines (specifically for lesion experimnets)  aLSO pUT IN VERTICAL LINES FOR important dates (e.g. lesions)
        if isfield(all_days_all_analysis.summary_of_experiment,'inter_lesion_days')==1;
            disp('adding line for baseline value');
            baseline_values_temp=[];
            for i=1:all_days_all_analysis.summary_of_experiment.inter_lesion_days(1);
                baseline_values_temp=[baseline_values_temp all_days_all_analysis.data{i}.all_songs.transitions{j}.(analysis_type)'];
            end
            baseline_mean=median(baseline_values_temp);
            baseline_std=std(baseline_values_temp);
            line(xlim, [baseline_mean baseline_mean],'Color','r','LineWidth',0.4); % mean lines
            line(xlim, [baseline_mean+baseline_std baseline_mean+baseline_std],'LineStyle','--','Color',[0.7 0.1 0.1],'LineWidth',0.3); % variance lines
            line(xlim, [baseline_mean-baseline_std baseline_mean-baseline_std],'LineStyle','--','Color', [0.7 0.1 0.1],'LineWidth',0.3);
            
            disp('adding lines for lesion days');
            % lesion 1
            for i=1:length(all_days_all_analysis.summary_of_experiment.inter_lesion_days);
                
                line([all_days_all_analysis.summary_of_experiment.inter_lesion_days(i)+0.1, ...
                    all_days_all_analysis.summary_of_experiment.inter_lesion_days(i)+0.1],ylim,'Color',[0.2 0.6 0.2],'LineWidth',1.8);
            end
        end
        
        % BASELINE AND IMPORTANT DATES FOR WN EXPERIMENTS
        if isfield(all_days_all_analysis.summary_of_experiment,'WN_baseline_days')==1;
            disp('adding line for baseline value');
            baseline_values_temp=[];
            for i=1:all_days_all_analysis.summary_of_experiment.WN_baseline_days;
                baseline_values_temp=[baseline_values_temp all_days_all_analysis.data{i}.all_songs.(transitions{j}).(analysis_type)'];
            end
            baseline_mean=median(baseline_values_temp);
            baseline_std=std(baseline_values_temp);
            line(xlim, [baseline_mean baseline_mean],'Color','r','LineWidth',0.4); % mean lines
            line(xlim, [baseline_mean+baseline_std baseline_mean+baseline_std],'LineStyle','--','Color',[0.7 0.1 0.1],'LineWidth',0.3); % variance lines
            line(xlim, [baseline_mean-baseline_std baseline_mean-baseline_std],'LineStyle','--','Color', [0.7 0.1 0.1],'LineWidth',0.3);
            
            disp('adding lines for WN on, driving--> consolid, and off days');
            line([all_days_all_analysis.summary_of_experiment.WN_baseline_days-0.5, ...
                all_days_all_analysis.summary_of_experiment.WN_baseline_days-0.5],ylim,'Color',[0.2 0.6 0.2],'LineWidth',1.8);% WN start
            
            line([all_days_all_analysis.summary_of_experiment.WN_baseline_days+all_days_all_analysis.summary_of_experiment.WN_driving_days-0.5, ...
                all_days_all_analysis.summary_of_experiment.WN_baseline_days+all_days_all_analysis.summary_of_experiment.WN_driving_days-0.5],ylim,'Color',[0.2 0.6 0.2],'LineWidth',1.8); % WN driving--> consolidation
            
            line([all_days_all_analysis.summary_of_experiment.WN_baseline_days+all_days_all_analysis.summary_of_experiment.WN_driving_days+all_days_all_analysis.summary_of_experiment.WN_consolidation_days-0.5, ...
                all_days_all_analysis.summary_of_experiment.WN_baseline_days+all_days_all_analysis.summary_of_experiment.WN_driving_days+all_days_all_analysis.summary_of_experiment.WN_consolidation_days-0.5],ylim,'Color',[0.2 0.6 0.2],'LineWidth',1.8);  % WN end
        end
        
        title([analysis_type ' over days for transition: ' transitions{j} '; dots and lines are medians']); xlabel('days');
    end
    
end

%% PLOT OVER ALL DAYS - for individual syllable analyses (plotting across days)
% can't remember if this is any analysis or just pitch.  if any, then how
% different from above.  if pitch, then how different from below
analysis_type='pitch';
syllables = 'fhi';
first_day_num=datenum(all_days_all_analysis.parameters.first_day);


for j=1:length(syllables);
    figure(j); hold on;
    for i=1:length(all_days_all_analysis.data); % number of days
        try
            scatter(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type)(:,1)-first_day_num,all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type)(:,2),'b');
        catch err
            continue
        end
    end
    
    % Put in medians, std
    for i=1:length(all_days_all_analysis.data); % number of days
        try
            errorbar(median(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type)(:,1)-first_day_num),median(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type)(:,2)),...
                std(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type)(:,2)),'Marker','.','MarkerSize',27,'Color','k');
        catch err
            continue
        end
    end
    
    % put in baseline lines (specifically for lesion experimnets, and PUT IN VERTICAL LINES FOR important dates (e.g. lesions)
    
    if isfield(all_days_all_analysis.summary_of_experiment,'inter_lesion_days')==1;
        disp('adding line for baseline value');
        baseline_values_temp=[];
        for i=1:all_days_all_analysis.summary_of_experiment.inter_lesion_days(1);
            baseline_values_temp=[baseline_values_temp all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type)(:,2)'];
        end
        baseline_mean=median(baseline_values_temp);
        baseline_std=std(baseline_values_temp);
        line(xlim, [baseline_mean baseline_mean],'Color','r','LineWidth',0.4); % mean lines
        line(xlim, [baseline_mean+baseline_std baseline_mean+baseline_std],'LineStyle','--','Color',[0.7 0.1 0.1],'LineWidth',0.3); % variance lines
        line(xlim, [baseline_mean-baseline_std baseline_mean-baseline_std],'LineStyle','--','Color', [0.7 0.1 0.1],'LineWidth',0.3);
        
        disp('adding lines for lesion days');
        % lesion 1
        for i=1:length(all_days_all_analysis.summary_of_experiment.inter_lesion_days);
            
            line([all_days_all_analysis.summary_of_experiment.inter_lesion_days(i)+0.1, ...
                all_days_all_analysis.summary_of_experiment.inter_lesion_days(i)+0.1],ylim,'Color',[0.2 0.6 0.2],'LineWidth',1.8);
        end
    end
    % MODIFY BELOW FROM TRANSITION TO SYLLABLE INPUTS
%     % BASELINE AND IMPORTANT DATES FOR WN EXPERIMENTS
%     if isfield(all_days_all_analysis.summary_of_experiment,'WN_baseline_days')==1;
%         disp('adding line for baseline value');
%         baseline_values_temp=[];
%         for i=1:all_days_all_analysis.summary_of_experiment.WN_baseline_days;
%             baseline_values_temp=[baseline_values_temp all_days_all_analysis.data{i}.all_songs.(transitions{j}).(analysis_type)'];
%         end
%         baseline_mean=median(baseline_values_temp);
%         baseline_std=std(baseline_values_temp);
%         line(xlim, [baseline_mean baseline_mean],'Color','r','LineWidth',0.4); % mean lines
%         line(xlim, [baseline_mean+baseline_std baseline_mean+baseline_std],'LineStyle','--','Color',[0.7 0.1 0.1],'LineWidth',0.3); % variance lines
%         line(xlim, [baseline_mean-baseline_std baseline_mean-baseline_std],'LineStyle','--','Color', [0.7 0.1 0.1],'LineWidth',0.3);
%         
%         disp('adding lines for WN on, driving--> consolid, and off days');
%         line([all_days_all_analysis.summary_of_experiment.WN_baseline_days-0.5, ...
%             all_days_all_analysis.summary_of_experiment.WN_baseline_days-0.5],ylim,'Color',[0.2 0.6 0.2],'LineWidth',1.8);% WN start
%         
%         line([all_days_all_analysis.summary_of_experiment.WN_baseline_days+all_days_all_analysis.summary_of_experiment.WN_driving_days-0.5, ...
%             all_days_all_analysis.summary_of_experiment.WN_baseline_days+all_days_all_analysis.summary_of_experiment.WN_driving_days-0.5],ylim,'Color',[0.2 0.6 0.2],'LineWidth',1.8); % WN driving--> consolidation
%         
%         line([all_days_all_analysis.summary_of_experiment.WN_baseline_days+all_days_all_analysis.summary_of_experiment.WN_driving_days+all_days_all_analysis.summary_of_experiment.WN_consolidation_days-0.5, ...
%             all_days_all_analysis.summary_of_experiment.WN_baseline_days+all_days_all_analysis.summary_of_experiment.WN_driving_days+all_days_all_analysis.summary_of_experiment.WN_consolidation_days-0.5],ylim,'Color',[0.2 0.6 0.2],'LineWidth',1.8);  % WN end
%     end
    
    title([analysis_type ' over days for syllable: ' syllables(j) '; dots and lines are medians']); xlabel('days');
end

%% below for PITCH (using db_contour)
% SEEMS SAME AS ABOVE - CAN'T REMEMBER HOW DIFFERENT...

if strcmp(analysis_type,'pitch')==1;
    
    for j=1:length(syllables);
        figure(j); hold on;
        % this is to make day 1 = day 0.
        
        NumToAddToAllDays= -mean(all_days_all_analysis.data{1}.all_songs.(syllables(j)).(analysis_type).FF(:,1));
        
        for i=1:length(all_days_all_analysis.data); % number of days
            try
                scatter(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type).FF(:,1)+NumToAddToAllDays,all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type).FF(:,2),'b');
            catch err
                continue
            end
        end
        
        % Put in medians, std
    for i=1:length(all_days_all_analysis.data); % number of days
        try
            errorbar(median(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type).FF(:,1))+NumToAddToAllDays,median(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type).FF(:,2)),...
                std(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type).FF(:,2)),'Marker','.','MarkerSize',27,'Color','k');
        catch err
            continue
        end
    end
    
    % put in baseline lines (specifically for lesion experimnets, and PUT IN VERTICAL LINES FOR important dates (e.g. lesions)
    
    if isfield(all_days_all_analysis.summary_of_experiment,'inter_lesion_days')==1;
        disp('adding line for baseline value');
        baseline_values_temp=[];
        for i=1:all_days_all_analysis.summary_of_experiment.inter_lesion_days(1);
            baseline_values_temp=[baseline_values_temp all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type).FF(:,2)'];
        end
        baseline_mean=median(baseline_values_temp);
        baseline_std=std(baseline_values_temp);
        line(xlim, [baseline_mean baseline_mean],'Color','r','LineWidth',0.4); % mean lines
        line(xlim, [baseline_mean+baseline_std baseline_mean+baseline_std],'LineStyle','--','Color',[0.7 0.1 0.1],'LineWidth',0.3); % variance lines
        line(xlim, [baseline_mean-baseline_std baseline_mean-baseline_std],'LineStyle','--','Color', [0.7 0.1 0.1],'LineWidth',0.3);
        
        disp('adding lines for lesion days');
        % lesion 1
        %         for i=1:length(all_days_all_analysis.summary_of_experiment.inter_lesion_days);
        %             line([all_days_all_analysis.summary_of_experiment.inter_lesion_days(i)+0.1, ...
        %                 all_days_all_analysis.summary_of_experiment.inter_lesion_days(i)+0.1],ylim,'Color',[0.2 0.6 0.2],'LineWidth',1.8);
        %         end
        for i=1:length(all_days_all_analysis.summary_of_experiment.inter_lesion_days); % USE THE GREEN VERSION ABOVE IF YOU KEEP AS DAYS IN BETWEEN LESIONS
            absolute_lesion_day= sum(all_days_all_analysis.summary_of_experiment.inter_lesion_days(1:i)); % i.e. 3 pre-lesion days implies lesion 1 on day 3 (since day 1 is really 0)
            line([absolute_lesion_day+0.1, ...
                absolute_lesion_day+0.1],ylim,'Color',[0.2 0.6 0.2],'LineWidth',1.8);
        end
    end
    % MODIFY BELOW FROM TRANSITION TO SYLLABLE INPUTS
    %     % BASELINE AND IMPORTANT DATES FOR WN EXPERIMENTS
    %     if isfield(all_days_all_analysis.summary_of_experiment,'WN_baseline_days')==1;
    %         disp('adding line for baseline value');
    %         baseline_values_temp=[];
    %         for i=1:all_days_all_analysis.summary_of_experiment.WN_baseline_days;
    %             baseline_values_temp=[baseline_values_temp all_days_all_analysis.data{i}.all_songs.(transitions{j}).(analysis_type)'];
%         end
%         baseline_mean=median(baseline_values_temp);
%         baseline_std=std(baseline_values_temp);
%         line(xlim, [baseline_mean baseline_mean],'Color','r','LineWidth',0.4); % mean lines
%         line(xlim, [baseline_mean+baseline_std baseline_mean+baseline_std],'LineStyle','--','Color',[0.7 0.1 0.1],'LineWidth',0.3); % variance lines
%         line(xlim, [baseline_mean-baseline_std baseline_mean-baseline_std],'LineStyle','--','Color', [0.7 0.1 0.1],'LineWidth',0.3);
%         
%         disp('adding lines for WN on, driving--> consolid, and off days');
%         line([all_days_all_analysis.summary_of_experiment.WN_baseline_days-0.5, ...
%             all_days_all_analysis.summary_of_experiment.WN_baseline_days-0.5],ylim,'Color',[0.2 0.6 0.2],'LineWidth',1.8);% WN start
%         
%         line([all_days_all_analysis.summary_of_experiment.WN_baseline_days+all_days_all_analysis.summary_of_experiment.WN_driving_days-0.5, ...
%             all_days_all_analysis.summary_of_experiment.WN_baseline_days+all_days_all_analysis.summary_of_experiment.WN_driving_days-0.5],ylim,'Color',[0.2 0.6 0.2],'LineWidth',1.8); % WN driving--> consolidation
%         
%         line([all_days_all_analysis.summary_of_experiment.WN_baseline_days+all_days_all_analysis.summary_of_experiment.WN_driving_days+all_days_all_analysis.summary_of_experiment.WN_consolidation_days-0.5, ...
%             all_days_all_analysis.summary_of_experiment.WN_baseline_days+all_days_all_analysis.summary_of_experiment.WN_driving_days+all_days_all_analysis.summary_of_experiment.WN_consolidation_days-0.5],ylim,'Color',[0.2 0.6 0.2],'LineWidth',1.8);  % WN end
%     end
    
    title([analysis_type ' over days for syllable: ' syllables(j) '; dots and lines are medians']); xlabel('days');
end
end
    

%% PITCH - analyzing pitch of a syllable in a specific context (e.g. b only in ab, not in db).
% Plot data from all_days_various.lt_db_check_templ_freq
% (see analysis for bk24bk63)
 
% PITCH COV (see analysis for bk24bk63)


%% annotating WN epochs and pitch thresholds - see bk24bk63 analysis.


%% to plot histogram of all vocalization durations.
day1=all_days_all_analysis.summary_of_experiment.inter_lesion_days(3);
day2=all_days_all_analysis.summary_of_experiment.inter_lesion_days(3)+1;
nbins=100;

figure;
hist(all_days_all_analysis.data{day1}.all_songs.all_vocalization_stats.durations_all_songs{day1},nbins);
figure;
hist(all_days_all_analysis.data{day2}.all_songs.all_vocalization_stats.durations_all_songs{day2},nbins);



%% To plot one variable against another (e.g. is my entropy increase real or does it just reflect changes in amplitude?)
analysis_type_one='amplitude'; % x axis
analysis_type_two='entropy'; % y axis
% syllables='cdabjklm';
syllables=all_days_all_analysis.parameters.syllables_seq_func;
syllables='fk';
lin_reg_on_or_off='off';
type_of_manipulation='WN';


switch type_of_manipulation;
    case 'WN'
        baseline_days=all_days_all_analysis.summary_of_experiment.WN_baseline_days;
        important_days=[all_days_all_analysis.summary_of_experiment.WN_baseline_days all_days_all_analysis.summary_of_experiment.WN_baseline_days+all_days_all_analysis.summary_of_experiment.WN_driving_days];
    case 'lesion'
        baseline_days=all_days_all_analysis.summary_of_experiment.inter_lesion_days(1);
        important_days=all_days_all_analysis.summary_of_experiment.inter_lesion_days;
end
    
    
% VERSION 1: look at all renditions for specific days (same plot)
for j=1:length(syllables);
    figure; hold on
    for i=1:baseline_days;
        try
            if strcmp(lin_reg_on_or_off,'on')==1;
                [r,m,b]=regression(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_one)(:,2), all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_two)(:,2),'one');
                plot(min(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_one)(:,2)):max(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_one)(:,2)),b+m*[min(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_one)(:,2)):max(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_one)(:,2))],'b');
            end
            scatter(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_one)(:,2), all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_two)(:,2),'b');
        catch err
            continue
        end
        
    end
    
    for i=important_days(2)-25:important_days(2)-15;
        try
            if strcmp(lin_reg_on_or_off,'on')==1;
                [r,m,b]=regression(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_one)(:,2), all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_two)(:,2),'one');
                plot(min(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_one)(:,2)):max(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_one)(:,2)),b+m*[min(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_one)(:,2)):max(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_one)(:,2))],'r');
            end
            scatter(all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_one)(:,2), all_days_all_analysis.data{i}.all_songs.(syllables(j)).(analysis_type_two)(:,2),'r');
        catch err
            continue
        end
    end
    title(['individual renditions of syllable ' syllables(j) '. ' analysis_type_two ' vs. ' analysis_type_one '. Pre=blue Post=red'])
    xlabel(analysis_type_one); ylabel(analysis_type_two);
end

% to save - incomplete
cd(['/bluejay2/lucas/birds/' all_days_all_analysis.parameters.name_of_bird '/all_days_all_analysis/PLOT'])

%----------------------------------------------------------------
% VERSION 2: to plot across days (one scatter point per day)

num_of_days=length(all_days_all_analysis.data);
divisions=1/num_of_days;

for ii=1:length(syllables);
    figure; hold on; title(['mean per day of syllable ' syllables(ii) '. ' analysis_type_two ' vs. ' analysis_type_one '. Pre=blue Post=green to red (chronological)'])
    for jj=1:num_of_days;
        try
            mean_per_day.(analysis_type_one){jj}=mean(all_days_all_analysis.data{jj}.all_songs.(syllables(ii)).(analysis_type_one)(:,2));
            mean_per_day.(analysis_type_two){jj}=mean(all_days_all_analysis.data{jj}.all_songs.(syllables(ii)).(analysis_type_two)(:,2));
            if jj<=baseline_days;
                plot(mean_per_day.(analysis_type_one){jj},mean_per_day.(analysis_type_two){jj},'.','MarkerSize',15,'Color','b')
            elseif jj>baseline_days;
                plot(mean_per_day.(analysis_type_one){jj},mean_per_day.(analysis_type_two){jj},'.','MarkerSize',15,'Color',jj*[divisions 0 0]+[0 1 0]-jj*[0 divisions 0])
            else
                disp('error')
            end
        catch
        end
    end
end

% FIND MEAN ACROSS DAYS (e.g. early, late learning)


%% analyzing data from all labels (e.g. motif durations) - see bk24bk63 analysis.


%% Plotting renditions/song for any syllable (used for rd66gr)

% change parameters:
syllables={'a','b','j'};

plot_colors=lt_make_plot_colors(size(syllables,2));

figure; hold on
for s=1:size(syllables,2);
    for j=1:length(all_days_all_analysis.data); % num of days
        
        syl=syllables{s};
        rend_per_song.(syl){j}=all_days_all_analysis.data{j}.summary_of_day.(syl).syl_rendition_amount...
            /all_days_all_analysis.data{j}.summary_of_day.song_amount_labeled;
        try
            hscatter3(s)= scatter(j,rend_per_song.(syl){j},'k','MarkerFaceColor',plot_colors{s});
        catch err
        end
    end
end

legend(hscatter3,syllables)
title('# renditions of each syl per song, across days'), xlabel('days'); ylabel('renditions/song');
