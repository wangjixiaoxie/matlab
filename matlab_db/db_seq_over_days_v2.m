% db_seq_over_days
%
% Written by DM Brady. 09/2012
%
% Takes .mat files produced by db_mult_seq_analy in the all_days_sequence
% folder for the bird and makes a series of timelines and calculations of
% learning

%% Collects user input on some basic information

parameters.birdname = input('What is the name of your bird?  ', 's');
parameters.phrase = input('Phrase?  ','s');
if ~isempty(strfind(parameters.phrase,'''')) == 1
    parameters.phrase =[];
end
parameters.con_or_div = input('Convergence or divergence? (con or div)  ','s');
parameters.day_begin = input('How many days before white noise?  ');
parameters.day_duration = input('How many days of white noise delivery?  ');


%% Loads list of .mat files

fid = fopen(['all_days_sequence_' parameters.phrase '/' parameters.con_or_div],'r');
current_line = fgetl(fid);
i=1;
while ischar(current_line)
    parameters.file_names(i,:) = current_line;
    parameters.dates(i,:) = current_line(length(parameters.birdname)+2:length(parameters.birdname)+10);
    current_line = fgetl(fid);
    i = i+1;
end
    

%% Timeline graphs

%loads multiple_trans from that day
for i = 1:size(parameters.file_names,1)
    load(['all_days_sequence_' parameters.phrase '/' parameters.file_names(i,:)])

    
    %% making the transitional probabilites timeline graph
    figure(1)
    hold on
    for j = 1:length(multiple_trans.parameters.tpofinterest)
        %plots transition probability
        h1(i,j) = plot(datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1,...
            multiple_trans.trans_prob.median.(multiple_trans.parameters.tpofinterest{j}),...
            'o',...
            'MarkerSize',10,...
            'Color', [1-((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))...
            0 0+((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))],...
            'MarkerFaceColor', [1-((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))...
            0 0+((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))]);
        
        %plots confidence interval
        line([datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1 datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1],...
            [multiple_trans.trans_prob.CI.(multiple_trans.parameters.tpofinterest{j})(1)...
            multiple_trans.trans_prob.CI.(multiple_trans.parameters.tpofinterest{j})(2)],...
            'LineWidth', 1, ...
            'Color', [1-((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))...
            0 0+((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))]);
        
        %plots line connecting to previous day's median
        current.trans_prob.trans_prob(i,j) = multiple_trans.trans_prob.median.(multiple_trans.parameters.tpofinterest{j});
        current.trans_prob.trans_prob(i,length(multiple_trans.parameters.tpofinterest)+1) = datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1;
        if i >= 2
            line([current.trans_prob.trans_prob(i-1,length(multiple_trans.parameters.tpofinterest)+1) current.trans_prob.trans_prob(i,length(multiple_trans.parameters.tpofinterest)+1)],...
                [current.trans_prob.trans_prob(i-1,j) current.trans_prob.trans_prob(i,j)],...
                'LineWidth', 2,...
                'Color', [1-((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))...
            0 0+((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))]);
        else
        end
    end
    
    %adds label, legend, and title to figure
    if i == size(parameters.file_names,1)
        title(['Transitional probabilities for ' parameters.birdname])
        xlabel('Day')
        ylabel('Probability')
        legend(h1(end,:), multiple_trans.parameters.tpofinterest)
        ylim([-0.1 1.1])
        xlim([0 datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+2])
        
        %line when white noise began
        line([parameters.day_begin+0.5 parameters.day_begin+0.5], [-0.1 1.1],...
            'LineStyle', '--',...
            'LineWidth', 1,...
            'Color', [0 0 0])
        
        %line when white noise ended
        line([parameters.day_begin+parameters.day_duration+0.5 parameters.day_begin+parameters.day_duration+0.5], [-0.1 1.1],...
            'LineStyle', '--',...
            'LineWidth', 1,...
            'Color', [0 0 0])
        
        %line for 0 and 1 probability
        line(get(gca,'XLim'),[0 0],...
            'LineStyle', '--',...
            'LineWidth', 1,...
            'Color', [.5 .5 .5])
        
        line(get(gca,'XLim'),[1 1],...
            'LineStyle', '--',...
            'LineWidth', 1,...
            'Color', [.5 .5 .5])
    else
    end
    hold off
    
    %% Standard deviation graph
    figure(2)
    hold on
    %Plots median standard deviation for that day
    plot(datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1,...
        multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{1}).sd,...
        'o',...
        'MarkerSize',10,...
        'Color', [0 0 0],...
        'MarkerFaceColor', [0 0 0]);
    
    %plots CI of standard deviation
    line([datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1 datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1],...
        [multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{1}).sd_CI(1)...
        multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{1}).sd_CI(2)],...
        'LineWidth', 1, ...
        'Color', [0 0 0]);
    
    %connects sd with previous day's sd
    current.sd.sd(i,1) = multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{1}).sd;
    current.sd.sd(i,2) = datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1;
    if i >= 2
            line([current.sd.sd(i-1,2) current.sd.sd(i,2)],...
                [current.sd.sd(i-1,1) current.sd.sd(i,1)],...
                'LineWidth', 2,...
                'Color', [0 0 0]);
    else
    end
    
    if i == size(parameters.file_names,1)
        title(['Variance in transitional probabilities across songs per day for ' parameters.birdname])
        xlabel('Day')
        ylabel('Standard Deviation')
        xlim([0 datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+2])
        h2 = get(gca,'Ylim');
        
        %line when white noise began
        line([parameters.day_begin+0.5 parameters.day_begin+0.5], [h2(1) h2(2)],...
            'LineStyle', '--',...
            'LineWidth', 1,...
            'Color', [0 0 0])
        
        %line when white noise ended
        line([parameters.day_begin+parameters.day_duration+0.5 parameters.day_begin+parameters.day_duration+0.5], [h2(1) h2(2)],...
            'LineStyle', '--',...
            'LineWidth', 1,...
            'Color', [0 0 0])
    else
    end
    hold off
    
    
    
%     %% Coefficient of variation
%     figure(3)
%     hold on
%     for j = 1:length(multiple_trans.parameters.tpofinterest)
%         %plots CV
%         h3(i,j) = plot(datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1,...
%             multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{j}).cv,...
%             'o',...
%             'MarkerSize',10,...
%             'Color', [1-((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))...
%             0 0+((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))],...
%             'MarkerFaceColor', [1-((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))...
%             0 0+((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))]);
%         
%         %plots confidence interval
%         line([datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1 datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1],...
%             [multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{j}).cv_CI(1)...
%             multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{j}).cv_CI(2)],...
%             'LineWidth', 1, ...
%             'Color', [1-((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))...
%             0 0+((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))]);
%         
%         %plots line connecting to previous day's median
%         current.CV.CV(i,j) = multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{j}).cv;
%         current.CV.CV(i,length(multiple_trans.parameters.tpofinterest)+1) = datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1;
%         if i >= 2
%             line([current.CV.CV(i-1,length(multiple_trans.parameters.tpofinterest)+1) current.CV.CV(i,length(multiple_trans.parameters.tpofinterest)+1)],...
%                 [current.CV.CV(i-1,j) current.CV.CV(i,j)],...
%                 'LineWidth', 2,...
%                 'Color', [1-((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))...
%             0 0+((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))]);
%         else
%         end
%     end
%     
%     %adds label, legend, and title to figure
%     if i == size(parameters.file_names,1)
%         title(['CV of transitional probabilities for ' parameters.birdname])
%         xlabel('Day')
%         ylabel('Coefficient of Variation')
%         legend(h3(end,:), multiple_trans.parameters.tpofinterest)
%         xlim([0 datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+2])
%         h4 = get(gca,'YLim');
%         
%         %line when white noise began
%         line([parameters.day_begin+0.5 parameters.day_begin+0.5], [h4(1) h4(2)],...
%             'LineStyle', '--',...
%             'LineWidth', 1,...
%             'Color', [0 0 0])
%         
%         %line when white noise ended
%         line([parameters.day_begin+parameters.day_duration+0.5 parameters.day_begin+parameters.day_duration+0.5], [h4(1) h4(2)],...
%             'LineStyle', '--',...
%             'LineWidth', 1,...
%             'Color', [0 0 0])
%     else
%     end
%     hold off
    
    
    %% Entropy
    figure(3)
    hold on
    %Plots median entropy for that day
    plot(datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1,...
        multiple_trans.entropy.median,...
        'o',...
        'MarkerSize',10,...
        'Color', [0 0 0],...
        'MarkerFaceColor', [0 0 0]);
    
    %plots CI of standard deviation
    line([datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1 datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1],...
        [multiple_trans.entropy.CI(1)...
        multiple_trans.entropy.CI(2)],...
        'LineWidth', 1, ...
        'Color', [0 0 0]);
    
    %connects sd with previous day's sd
    current.ent.ent(i,1) = multiple_trans.entropy.median;
    current.ent.ent(i,2) = datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1;
    if i >= 2
            line([current.ent.ent(i-1,2) current.ent.ent(i,2)],...
                [current.ent.ent(i-1,1) current.ent.ent(i,1)],...
                'LineWidth', 2,...
                'Color', [0 0 0]);
    else
    end
    
    if i == size(parameters.file_names,1)
        title(['Entropy per day for ' parameters.birdname])
        xlabel('Day')
        ylabel('Entropy (bits)')
        xlim([0 datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+2])
        h5 = get(gca,'Ylim');
        
        %line when white noise began
        line([parameters.day_begin+0.5 parameters.day_begin+0.5], [h5(1) h5(2)],...
            'LineStyle', '--',...
            'LineWidth', 1,...
            'Color', [0 0 0])
        
        %line when white noise ended
        line([parameters.day_begin+parameters.day_duration+0.5 parameters.day_begin+parameters.day_duration+0.5], [h5(1) h5(2)],...
            'LineStyle', '--',...
            'LineWidth', 1,...
            'Color', [0 0 0])
    else
    end
    hold off
    
    %% History dependence
    figure(4)
    hold on
    for j = 1:length(multiple_trans.parameters.tpofinterest)
        subplot(length(multiple_trans.parameters.tpofinterest),1,j)
        hold on
        if multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{j})(1) >...
                multiple_trans.within_shuffled.CI.(multiple_trans.parameters.tpofinterest{j})(2)
            current.hist_dep.color{i}{j} = [1 0 0];
        elseif multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{j})(1) == 1 && ...
                multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{j})(1) == 2
            current.hist_dep.color{i}{j} = [1 0 0];
        else
            current.hist_dep.color{i}{j} = [0 0 0];
        end
        
        %plots history dependence for a transition
        h6 = plot(datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1,...
            multiple_trans.hist_dep.median.(multiple_trans.parameters.tpofinterest{j}),...
            '<',...
            'MarkerSize',10,...
            'Color', current.hist_dep.color{i}{j},...
            'MarkerFaceColor', current.hist_dep.color{i}{j});
        
        %plots confidence interval of history dependence
        line([datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1 datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1],...
            [multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{j})(1)...
            multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{j})(2)],...
            'LineWidth', 2, ...
            'Color', current.hist_dep.color{i}{j});
        
        try
            %plots confidence interval of shuffled history dependence
            h7 = line([datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+0.75 datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+0.75],...
                [multiple_trans.within_shuffled.CI.(multiple_trans.parameters.tpofinterest{j})(1)...
                multiple_trans.within_shuffled.CI.(multiple_trans.parameters.tpofinterest{j})(2)],...
                'LineWidth', 2, ...
                'Color', [0.6 0.6 0.6]);
            
            if i == size(parameters.file_names,1)
                %label of subplot
                title(['History dependence of ' multiple_trans.parameters.tpofinterest{j}])
                ylabel('History dependence |p(ab|ab) - p(ab|ac)|')
                xlabel('Day')
                %             legend(h6, multiple_trans.parameters.tpofinterest{j})
                %             legend(h7, 'shuffled CI')
                h8 = get(gca,'Ylim');
                
                %line when white noise began
                line([parameters.day_begin+0.5 parameters.day_begin+0.5], [h8(1) h8(2)],...
                    'LineStyle', '--',...
                    'LineWidth', 1,...
                    'Color', [0 0 0])
                
                %line when white noise ended
                line([parameters.day_begin+parameters.day_duration+0.5 parameters.day_begin+parameters.day_duration+0.5], [h8(1) h8(2)],...
                    'LineStyle', '--',...
                    'LineWidth', 1,...
                    'Color', [0 0 0])
            else
            end
        catch err
            continue
        end
        
        
%         %plots line connecting to previous day's median
        current.hist_dep.hist_dep(i,j) = multiple_trans.hist_dep.median.(multiple_trans.parameters.tpofinterest{j});
        current.hist_dep.shuffled{i}{j} = multiple_trans.within_shuffled.CI.(multiple_trans.parameters.tpofinterest{j});
        current.hist_dep.hist_dep(i,length(multiple_trans.parameters.tpofinterest)+1) = datenum(parameters.dates(i,:))-datenum(parameters.dates(1,:))+1;
%         if i >= 2
%             line([current.trans_prob(i-1,length(multiple_trans.parameters.tpofinterest)+1) current.trans_prob(i,length(multiple_trans.parameters.tpofinterest)+1)],...
%                 [current.trans_prob(i-1,j) current.trans_prob(i,j)],...
%                 'LineWidth', 2,...
%                 'Color', [1-((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))...
%             0 0+((length(multiple_trans.parameters.tpofinterest)-j)./length(multiple_trans.parameters.tpofinterest))]);
%         else
%         end
    end
        
        
    
    clear multiple_trans
end

%gets rids of figure handles and local variables
clear('h*','i','j')

 

%% Calculations over WN experiment

%% Pulls out relevant values from 'current' structure

parameters.fieldnames_current = fieldnames(current);
for i = 1:length(fieldnames(current))
        calculations.(parameters.fieldnames_current{i}).baseline =...
            nanmean(current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(1:parameters.day_begin,:),1);
        
        calculations.(parameters.fieldnames_current{i}).WN =...
            nanmean(current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(parameters.day_begin:parameters.day_begin+parameters.day_duration,:),1);
        
        calculations.(parameters.fieldnames_current{i}).post =...
            nanmean(current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(parameters.day_begin+parameters.day_duration+1:end,:),1);
        
        calculations.(parameters.fieldnames_current{i}).max_WN =...
            max(current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(parameters.day_begin:parameters.day_begin+parameters.day_duration,:));
        calculations.(parameters.fieldnames_current{i}).min_WN =...
            min(current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(parameters.day_begin:parameters.day_begin+parameters.day_duration,:));
        
        calculations.(parameters.fieldnames_current{i}).final_WN =...
            current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(parameters.day_begin+parameters.day_duration,:);
        
        calculations.(parameters.fieldnames_current{i}).max_postWN = ...
            max(current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(parameters.day_begin+parameters.day_duration:end,:));
        calculations.(parameters.fieldnames_current{i}).min_postWN =...
            min(current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(parameters.day_begin+parameters.day_duration:end,:));
        
        try
            calculations.(parameters.fieldnames_current{i}).post_3d = ...
                current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(parameters.day_begin+parameters.day_duration+3,:);
            
            calculations.(parameters.fieldnames_current{i}).post_7d = ...
                current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(parameters.day_begin+parameters.day_duration+7,:);
        catch err
            continue
        end

end

%% Drawing baseline lines in figures

for i = 1:length(parameters.fieldnames_current)
    if strcmpi(parameters.fieldnames_current{i},'hist_dep') ~= 1
        figure(i)
        hold on
        h = get(gca,'XLim');
        for j = 1:size(calculations.(parameters.fieldnames_current{i}).baseline,2)-1
            line([h(1) h(2)], [calculations.(parameters.fieldnames_current{i}).baseline(j)...
                calculations.(parameters.fieldnames_current{i}).baseline(j)],...
                'LineStyle','--',...
                'Color', [.5 .5 .5]);
        end
    elseif strcmpi(parameters.fieldnames_current{i},'hist_dep') == 1
        figure(i)
        hold on
        for j = 1:length(findall(gcf,'type','axes'))
            subplot(length(findall(gcf,'type','axes')),1,j)
            hold on
            h = get(gca,'XLim');
            line([h(1) h(2)], [calculations.(parameters.fieldnames_current{i}).baseline(j)...
                calculations.(parameters.fieldnames_current{i}).baseline(j)],...
                'LineStyle','--',...
                'Color', [.5 .5 .5]);
        end
    end
end
        
clear('h*','j','i')            

%% Calculating percent changes from baseline (always from dominant transition)
        
parameters.dom_trans_col = find(max(calculations.trans_prob.baseline(:,1:end-1)) == calculations.trans_prob.baseline(:,1:end-1));

% Percent change between baseline and nth day 
for i = 1:length(fieldnames(current))
    calculations.(parameters.fieldnames_current{i}).percent_from_baseline =...
        (current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(:,parameters.dom_trans_col)-...
        calculations.(parameters.fieldnames_current{i}).baseline(:,parameters.dom_trans_col))./...
        calculations.(parameters.fieldnames_current{i}).baseline(:,parameters.dom_trans_col);
    if strcmpi(parameters.fieldnames_current{i},'trans_prob') == 1
        calculations.(parameters.fieldnames_current{i}).percent_from_baseline = -100*calculations.(parameters.fieldnames_current{i}).percent_from_baseline;
    else
        calculations.(parameters.fieldnames_current{i}).percent_from_baseline = 100*calculations.(parameters.fieldnames_current{i}).percent_from_baseline;
    end
end

% Correlation coefficient matrix between variables
for i = 1:length(parameters.fieldnames_current)
    calculations.correlation.matrix_baseline(:,i) = calculations.(parameters.fieldnames_current{i}).percent_from_baseline;
end

%Calculates correlation coefficient, p value, and 95% CI between different
%parameters
[calculations.correlation.baseline.r, calculations.correlation.baseline.p,...
    calculations.correlation.baseline.rlo, calculations.correlation.baseline.rup] = ...
    corrcoef(calculations.correlation.matrix_baseline);

calculations.correlation.baseline.names = parameters.fieldnames_current;

%Figures of trans prob and ent vs. history dependence
corr_interest.fieldnames = {'trans_prob' 'ent'};
corr_interest.labels = {'transition probability' 'entropy'};
for i = 1:length(corr_interest.fieldnames)
    figure(i+4)
    hold on
    %plot before learning
    h.before = plot(calculations.(corr_interest.fieldnames{i}).percent_from_baseline(1:parameters.day_begin),...
        calculations.hist_dep.percent_from_baseline(1:parameters.day_begin),...
        'ko',...
        'MarkerSize', 8,...
        'MarkerFaceColor', 'k');
    h.during = plot(calculations.(corr_interest.fieldnames{i}).percent_from_baseline(parameters.day_begin+1:parameters.day_begin+parameters.day_duration),...
        calculations.hist_dep.percent_from_baseline(parameters.day_begin+1:parameters.day_begin+parameters.day_duration),...
        'ro',...
        'MarkerSize', 8,...
        'MarkerFaceColor', 'r');
    h.after = plot(calculations.(corr_interest.fieldnames{i}).percent_from_baseline(parameters.day_begin+parameters.day_duration+1:end),...
        calculations.hist_dep.percent_from_baseline(parameters.day_begin+parameters.day_duration+1:end),...
        'bo',...
        'MarkerSize', 8,...
        'MarkerFaceColor', 'b');
    legend('before WN','during WN', 'after WN', 'Location', 'Best')
    title(['Change in ' corr_interest.labels{i} ' vs. change in history dependence from baseline'])
    if strcmpi(corr_interest.fieldnames{i}, 'trans_prob') == 1
        xlabel(['% Change in ' corr_interest.labels{i} ' from baseline (dominant transistion)'])
    else
        xlabel(['% Change in ' corr_interest.labels{i} ' from baseline'])
    end
    ylabel('% Change in History Dependence from baseline')
    h.y = get(gca,'Ylim');
    h.x = get(gca,'XLim');
    text(min(h.x)*(3/4),max(h.y)*(3/4),['r = ' ...
        num2str(calculations.correlation.baseline.r(strcmpi(calculations.correlation.baseline.names,corr_interest.fieldnames{i}),...
        end))]);
    text(min(h.x)*(3/4),max(h.y)*(2/3),['p = ' ...
        num2str(calculations.correlation.baseline.p(strcmpi(calculations.correlation.baseline.names,corr_interest.fieldnames{i}),...
        end))]);
end
clear('corr_interest','h')


%% Calculating percent changes from previous day (always from dominant transition)

% Percent change between current day and previous.
for i = 1:length(parameters.fieldnames_current)
    for j = 1:length(current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(:,parameters.dom_trans_col))-1
        calculations.(parameters.fieldnames_current{i}).percent_from_previous_day(j) = ...
            (current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(j+1,parameters.dom_trans_col)-...
            current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(j,parameters.dom_trans_col))./...
            current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(j,parameters.dom_trans_col);
    end
    if strcmpi(parameters.fieldnames_current{i},'trans_prob') == 1
        calculations.(parameters.fieldnames_current{i}).percent_from_previous_day = -100*calculations.(parameters.fieldnames_current{i}).percent_from_previous_day';
    else
        calculations.(parameters.fieldnames_current{i}).percent_from_previous_day = 100*calculations.(parameters.fieldnames_current{i}).percent_from_previous_day';
    end
end

% Correlation coefficient matrix between variables
for i = 1:length(parameters.fieldnames_current)
    calculations.correlation.matrix_previous(:,i) = calculations.(parameters.fieldnames_current{i}).percent_from_previous_day;
end

%Calculates correlation coefficient, p value, and 95% CI between different
%parameters
[calculations.correlation.previous.r, calculations.correlation.previous.p, ...
    calculations.correlation.previous.rlo, calculations.correlation.previous.rup] = ...
    corrcoef(calculations.correlation.matrix_previous);

calculations.correlation.previous.names = parameters.fieldnames_current;

% %% Calculations max %change and %recovery
% 
% % Percent change from baseline (max during WN) and day it happened
% for i = 1:length(parameters.fieldnames_current)
%     calculations.(parameters.fieldnames_current{i}).percent_max =...
%         [max(calculations.(parameters.fieldnames_current{i}).percent_from_baseline(parameters.day_begin:parameters.day_duration))... 
%         current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(...
%         min(find(find(max(calculations.(parameters.fieldnames_current{i}).percent_from_baseline(parameters.day_begin:parameters.day_duration))...
%         == calculations.(parameters.fieldnames_current{i}).percent_from_baseline) > parameters.day_begin)),end)];
% end
% 
% % Percent change from baseline on the last day of WN
% for i = 1:length(parameters.fieldnames_current)
%     calculations.(parameters.fieldnames_current{i}).percent_lastWN =...
%         [calculations.(parameters.fieldnames_current{i}).percent_from_baseline(parameters.day_begin+parameters.day_duration)... 
%         parameters.day_begin+parameters.day_duration];
% end
% 
% try
%     % Percent recovery on third and seventh day after white noise from max
%     for i = 1:length(parameters.fieldnames_current)
%         days = [3 7];
%         for j = 1:length(days)
%             calculations.(parameters.fieldnames_current{i}).(['percent_max_recovery_' num2str(days(j)) 'd']) = ...
%                 (current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(parameters.day_begin+parameters.day_duration+days(j),...
%                 parameters.dom_trans_col)-current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})...
%                 (calculations.(parameters.fieldnames_current{i}).percent_max(2),parameters.dom_trans_col))./...
%                 (calculations.(parameters.fieldnames_current{i}).baseline(parameters.dom_trans_col)-...
%                 current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})...
%                 (calculations.(parameters.fieldnames_current{i}).percent_max(2),parameters.dom_trans_col))*100;
%         end
%         clear days
%     end
%     
%     % Percent recovery on third and seventh day after white noise from max
%     for i = 1:length(parameters.fieldnames_current)
%         days = [3 7];
%         for j = 1:length(days)
%             calculations.(parameters.fieldnames_current{i}).(['percent_lastWN_recovery_' num2str(days(j)) 'd']) = ...
%                 (current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})(parameters.day_begin+parameters.day_duration+days(j),...
%                 parameters.dom_trans_col)-current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})...
%                 (parameters.day_begin+parameters.day_duration,parameters.dom_trans_col))./...
%                 (calculations.(parameters.fieldnames_current{i}).baseline(parameters.dom_trans_col)-...
%                 current.(parameters.fieldnames_current{i}).(parameters.fieldnames_current{i})...
%                 (parameters.day_begin+parameters.day_duration,parameters.dom_trans_col))*100;
%         end
%         clear days
%     end
% catch err
%     continue
% end
% 
% 
% %% Cumulative changes in learning 
% 
% % Cumulative to max and last day of WN
% for i = 1:length(parameters.fieldnames_current)
%     calculations.(parameters.fieldnames_current{i}).cumulative_maxWN = calculations.(parameters.fieldnames_current{i}).percent_from_baseline./...
%         max(calculations.(parameters.fieldnames_current{i}).percent_from_baseline)*100;
% 
%     calculations.(parameters.fieldnames_current{i}).cumulative_lastWN = calculations.(parameters.fieldnames_current{i}).percent_from_baseline./...
%         calculations.(parameters.fieldnames_current{i}).percent_from_baseline(parameters.day_begin+parameters.day_duration)*100;
% end
% 
% clear('i','j')
%         
%         
 %% Save figures and workspace variables

if exist(['all_days_sequence_' parameters.phrase '/seq_over_days'],'dir') ~= 7
    mkdir(['all_days_sequence_' parameters.phrase '/seq_over_days'])
end
    
 save(['all_days_sequence_' parameters.phrase '/seq_over_days' '/' parameters.birdname '_' parameters.con_or_div])
 
 saveas(figure(1), ['all_days_sequence_' parameters.phrase '/seq_over_days' '/Transition_probabilities_for_'...
     parameters.birdname '_' parameters.con_or_div '_over_time'], 'fig')
 saveas(figure(2), ['all_days_sequence_' parameters.phrase '/seq_over_days' '/SD_of_transition_probabilities_for_'...
     parameters.birdname '_' parameters.con_or_div '_over_time'], 'fig')
 % saveas(figure(3), [parameters.computer parameters.birdname...
 %     '/all_days_sequence/CV_of_Transition_probabilities_for_' parameters.birdname '_' parameters.con_or_div '_over_time'], 'fig')
 saveas(figure(3), ['all_days_sequence_' parameters.phrase '/seq_over_days' '/Entropy_for_'...
     parameters.birdname '_' parameters.con_or_div '_over_time'], 'fig')
 saveas(figure(4), ['all_days_sequence_' parameters.phrase '/seq_over_days' '/History_depedence_for_'...
     parameters.birdname '_' parameters.con_or_div '_over_time'], 'fig')
 saveas(figure(5), ['all_days_sequence_' parameters.phrase '/seq_over_days' '/Trans_prob_vs_hist_dep_for_'...
     parameters.birdname '_' parameters.con_or_div '_over_time'], 'fig')
 saveas(figure(6), ['all_days_sequence_' parameters.phrase '/seq_over_days' '/Ent_vs_hist_dep_for_'...
     parameters.birdname '_' parameters.con_or_div '_over_time'], 'fig')

        
        
        
        
        
        
    