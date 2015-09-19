% 
% %% make structures to know the timing and notes associated with triggers
% for i = 1:length(syllables.(tpofinterest{1}).match.trigger)
%     history_seq.(tpofinterest{1}){i} =...
%         syllables.(tpofinterest{1}).match.trigger(i).ttimes(syllables.(tpofinterest{1}).match.trigger(i).tnotes == tpofinterest{1}(2));
%     history_seq.(tpofinterest{2}){i} =...
%         syllables.(tpofinterest{2}).match.trigger(i).ttimes(syllables.(tpofinterest{2}).match.trigger(i).tnotes == tpofinterest{2}(2));
% end
% 
% %% combines trigger times for syllables of interest into one group and then
% %replaces them with the notes. Final product is a list per song of
% %syllables in order
% for i = 1:length(syllables.(tpofinterest{1}).match.trigger)
%     history_seq.time_order{i} = sort([history_seq.(tpofinterest{1}){i}; history_seq.(tpofinterest{2}){i}]);
%     for j = 1:length(history_seq.time_order{i})
%         if sum(history_seq.time_order{i}(j) == syllables.(tpofinterest{1}).match.trigger(i).ttimes) == 1
%             history_seq.syl_order{i}(j) = tpofinterest{1}(2);
%         elseif sum(history_seq.time_order{i}(j) == syllables.(tpofinterest{2}).match.trigger(i).ttimes) == 1
%             history_seq.syl_order{i}(j) = tpofinterest{2}(2);
%         end
%     end
% end

%% gets labels of songs
history_seq.syl_order = getlabels(batchfile);

%% gets rid unlabeled trials

for i  = 1:length(history_seq.syl_order)
    history_seq.syl_order{i} = history_seq.syl_order{i}(history_seq.syl_order{i} ~= '-');
    history_seq.syl_order{i} = history_seq.syl_order{i}(history_seq.syl_order{i} == tpofinterest{1}(end) |...
        history_seq.syl_order{i} == tpofinterest{2}(end));
end

history_seq.syl_order = history_seq.syl_order(~cellfun('isempty',history_seq.syl_order));

%% calculates transition probabilities
for i = 1:length(history_seq.syl_order)
    history_seq.transition_prob.([tpofinterest{1} '_to_' tpofinterest{1}]){i} = 0;
    history_seq.transition_prob.([tpofinterest{1} '_to_' tpofinterest{2}]){i} = 0;
    history_seq.transition_prob.([tpofinterest{2} '_to_' tpofinterest{1}]){i} = 0;
    history_seq.transition_prob.([tpofinterest{2} '_to_' tpofinterest{2}]){i} = 0;
    
    for j = 1:length(history_seq.syl_order{i})-1
        
        %for transition from ab to ab
        if history_seq.syl_order{i}(j) == tpofinterest{1}(2) && history_seq.syl_order{i}(j) == history_seq.syl_order{i}(j+1)
            history_seq.transition_prob.([tpofinterest{1} '_to_' tpofinterest{1}]){i} = ...
                history_seq.transition_prob.([tpofinterest{1} '_to_' tpofinterest{1}]){i} + 1;
            history_seq.transition_prob.([tpofinterest{1} '_to_' tpofinterest{2}]){i} = ...
                history_seq.transition_prob.([tpofinterest{1} '_to_' tpofinterest{2}]){i} + 0; 
            
        %for transition from ab to ac
        elseif history_seq.syl_order{i}(j) == tpofinterest{1}(2) && history_seq.syl_order{i}(j) ~= history_seq.syl_order{i}(j+1)
            history_seq.transition_prob.([tpofinterest{1} '_to_' tpofinterest{1}]){i} = ...
                history_seq.transition_prob.([tpofinterest{1} '_to_' tpofinterest{1}]){i} + 0;
            history_seq.transition_prob.([tpofinterest{1} '_to_' tpofinterest{2}]){i} = ...
                history_seq.transition_prob.([tpofinterest{1} '_to_' tpofinterest{2}]){i} + 1;
            
        %for transition from ac to ab
        elseif history_seq.syl_order{i}(j) == tpofinterest{2}(2) && history_seq.syl_order{i}(j) ~= history_seq.syl_order{i}(j+1)
            history_seq.transition_prob.([tpofinterest{2} '_to_' tpofinterest{1}]){i} = ...
                history_seq.transition_prob.([tpofinterest{2} '_to_' tpofinterest{1}]){i} + 1;
            history_seq.transition_prob.([tpofinterest{2} '_to_' tpofinterest{2}]){i} = ...
                history_seq.transition_prob.([tpofinterest{2} '_to_' tpofinterest{2}]){i} + 0;
        
        %for transition from ac to ac    
        elseif history_seq.syl_order{i}(j) == tpofinterest{2}(2) && history_seq.syl_order{i}(j) == history_seq.syl_order{i}(j+1)
            history_seq.transition_prob.([tpofinterest{2} '_to_' tpofinterest{1}]){i} = ...
               history_seq.transition_prob.([tpofinterest{2} '_to_' tpofinterest{1}]){i} + 0;
            history_seq.transition_prob.([tpofinterest{2} '_to_' tpofinterest{2}]){i} = ...
                history_seq.transition_prob.([tpofinterest{2} '_to_' tpofinterest{2}]){i} + 1;  
        end
    end
end

%creates a variable with the list of possible transitions.
trans_hist_interest = fieldnames(history_seq.transition_prob);

%calculates transition probability
for i = 1:length(trans_hist_interest)
    for j = 1:length(history_seq.transition_prob.(trans_hist_interest{i}))
        if mod(i,2) == 1
            history_seq.transition_prob.(trans_hist_interest{i}){j} = history_seq.transition_prob.(trans_hist_interest{i}){j}./...
                (history_seq.transition_prob.(trans_hist_interest{i}){j} + history_seq.transition_prob.(trans_hist_interest{i+1}){j});
        elseif mod(i,2) == 0
            history_seq.transition_prob.(trans_hist_interest{i}){j} = 1-history_seq.transition_prob.(trans_hist_interest{i-1}){j};
        end
    end
end

%gets rid of empty transition probability cells
    
for i = 1:length(trans_hist_interest)
    history_seq.transition_prob.(trans_hist_interest{i}) = ...
        history_seq.transition_prob.(trans_hist_interest{i})(~cellfun(@isnan,...
        history_seq.transition_prob.(trans_hist_interest{i})));
end

%converts transition probability structure from cell to matrix

for i = 1:length(trans_hist_interest)
    history_seq.transition_prob.(trans_hist_interest{i}) = ...
        cell2mat(history_seq.transition_prob.(trans_hist_interest{i}));
end


% bootstrap mean and 95% CI for history dependence

numbootstrptrials = 10000;

for i = 1:numbootstrptrials
    for j = 1:length(trans_hist_interest)
        
        history_seq.bootstrp.(trans_hist_interest{j}).sampling{i} = randi(size(history_seq.transition_prob.(trans_hist_interest{j}),2),...
            [1 size(history_seq.transition_prob.(trans_hist_interest{j}),2)]);
        
        history_seq.bootstrp.(trans_hist_interest{j}).boot_run.run{i} =...
            history_seq.transition_prob.(trans_hist_interest{j})(history_seq.bootstrp.(trans_hist_interest{j}).sampling{i});
        
        history_seq.bootstrp.(trans_hist_interest{j}).boot_run.mean(i) = nanmean(history_seq.bootstrp.(trans_hist_interest{j}).boot_run.run{i});
        
    end
end

%calculates median and 95% CI for conditional probabilities and makes
%histograms of bootstrap runs
for j = 1:length(trans_hist_interest)
    history_seq.bootstrp.(trans_hist_interest{j}).boot_run.median_of_mean = median(history_seq.bootstrp.(trans_hist_interest{j}).boot_run.mean);
    history_seq.bootstrp.(trans_hist_interest{j}).boot_run.CI_of_mean = ...
        [prctile(history_seq.bootstrp.(trans_hist_interest{j}).boot_run.mean,5)...
        prctile(history_seq.bootstrp.(trans_hist_interest{j}).boot_run.mean,95)];
    
%     %histograms of mean transition probabilities
%     [history_seq.bootstrp.(trans_hist_interest{j}).boot_run.hist.count history_seq.bootstrp.(trans_hist_interest{j}).boot_run.hist.bins] = ...
%         hist(history_seq.bootstrp.(trans_hist_interest{j}).boot_run.mean,100);
%     history_seq.bootstrp.(trans_hist_interest{j}).boot_run.hist.rel_freq = ...
%         history_seq.bootstrp.(trans_hist_interest{j}).boot_run.hist.count./sum(history_seq.bootstrp.(trans_hist_interest{j}).boot_run.hist.count);
%     
%     figure(12345+j)
%     hold on
%     bar(history_seq.bootstrp.(trans_hist_interest{j}).boot_run.hist.bins, history_seq.bootstrp.(trans_hist_interest{j}).boot_run.hist.rel_freq)
%     db_prepare_bootstrp_hist(...
%         ['Conditional probability for ' trans_hist_interest{j} ' on ' date_of_exp '  (bootstrap' num2str(numbootstrptrials) ' trials)'],...
%         'Relative Frequency',...
%         'Conditional Probability',...
%         history_seq.bootstrp.(trans_hist_interest{j}).boot_run.CI_of_mean,...
%         history_seq.bootstrp.(trans_hist_interest{j}).boot_run.hist.rel_freq,...
%         history_seq.bootstrp.(trans_hist_interest{j}).boot_run.median_of_mean)
%     
%     saveas(figure(12345+j), ['check_sequence_timing_' today_date '/Conditional_probability_' trans_hist_interest{j} '_' date_of_exp], 'fig')
end

%% calculating history depedence (ex: |p(ab|ab)-p(ab|ac)|)

for i = 1:length(tpofinterest)
    switch i
        case 1
            history_seq.hist_dep.(tpofinterest{i}).bootstrap.boot_run = abs(history_seq.bootstrp.(trans_hist_interest{1}).boot_run.mean-...
                history_seq.bootstrp.(trans_hist_interest{3}).boot_run.mean);
        case 2
            history_seq.hist_dep.(tpofinterest{i}).bootstrap.boot_run = abs(history_seq.bootstrp.(trans_hist_interest{2}).boot_run.mean-...
                history_seq.bootstrp.(trans_hist_interest{4}).boot_run.mean);
    end
    
    history_seq.hist_dep.(tpofinterest{i}).median = nanmedian(history_seq.hist_dep.(tpofinterest{i}).bootstrap.boot_run);
    history_seq.hist_dep.(tpofinterest{i}).CI = [prctile(history_seq.hist_dep.(tpofinterest{i}).bootstrap.boot_run,5)...
        prctile(history_seq.hist_dep.(tpofinterest{i}).bootstrap.boot_run,95)];
    
    %making histograms
    [history_seq.hist_dep.(tpofinterest{i}).bootstrap.hist.count history_seq.hist_dep.(tpofinterest{i}).bootstrap.hist.bins] = ...
        hist(history_seq.hist_dep.(tpofinterest{i}).bootstrap.boot_run,100);
    history_seq.hist_dep.(tpofinterest{i}).bootstrap.hist.rel_freq = ...
        history_seq.hist_dep.(tpofinterest{i}).bootstrap.hist.count./sum(history_seq.hist_dep.(tpofinterest{i}).bootstrap.hist.count);
    
%     figure(123456+i)
%     hold on
%     bar(history_seq.hist_dep.(tpofinterest{i}).bootstrap.hist.bins, history_seq.hist_dep.(tpofinterest{i}).bootstrap.hist.rel_freq)
%     db_prepare_bootstrp_hist(...
%         ['History depedence for ' tpofinterest{i} ' on ' date_of_exp '  (bootstrap' num2str(numbootstrptrials) ' trials)'],...
%         'Relative Frequency',...
%         'Conditional Probability',...
%         history_seq.hist_dep.(tpofinterest{i}).CI,...
%         history_seq.hist_dep.(tpofinterest{i}).bootstrap.hist.rel_freq,...
%         history_seq.hist_dep.(tpofinterest{i}).median)
    
    % saveas(figure(123456+i), ['check_sequence_timing_' today_date '/History_dependence_' tpofinterest{i} '_' date_of_exp], 'fig')
end

%saves history dependence bootstrap for only first branch point (they
%should be the same)
figure(123456)
    hold on
    bar(history_seq.hist_dep.(tpofinterest{1}).bootstrap.hist.bins, history_seq.hist_dep.(tpofinterest{1}).bootstrap.hist.rel_freq)
    db_prepare_bootstrp_hist(...
        ['History depedence for ' bird_name ' on ' date_of_exp '  (bootstrap' num2str(numbootstrptrials) ' trials)'],...
        'Relative Frequency',...
        'Conditional Probability',...
        history_seq.hist_dep.(tpofinterest{1}).CI,...
        history_seq.hist_dep.(tpofinterest{1}).bootstrap.hist.rel_freq,...
        history_seq.hist_dep.(tpofinterest{1}).median)
saveas(figure(123456), ['check_sequence_timing_' today_date '/History_dependence_' bird_name '_' date_of_exp], 'fig')


%% calculating shuffled songs (baseline)

%Shuffles songs
for j = 1:length(history_seq.syl_order)
    history_seq.shuffled_songs.syl_order{j} = cell(1,numbootstrptrials);
    for i = 1:numbootstrptrials
        history_seq.shuffled_songs.syl_order{j}{i} = history_seq.syl_order{j}(randperm(length(history_seq.syl_order{j})));
    end
end

%First step to calculating transistion probabilities

for i = 1:length(history_seq.shuffled_songs.syl_order)
    history_seq.shuffled_songs.(trans_hist_interest{1}){i} = cell(numbootstrptrials,1);
    history_seq.shuffled_songs.(trans_hist_interest{2}){i} = cell(numbootstrptrials,1);
    history_seq.shuffled_songs.(trans_hist_interest{3}){i} = cell(numbootstrptrials,1);
    history_seq.shuffled_songs.(trans_hist_interest{4}){i} = cell(numbootstrptrials,1);
    
    for j = 1:length(history_seq.shuffled_songs.syl_order{i})
        history_seq.shuffled_songs.(trans_hist_interest{1}){i}{j} = 0;
        history_seq.shuffled_songs.(trans_hist_interest{2}){i}{j} = 0;
        history_seq.shuffled_songs.(trans_hist_interest{3}){i}{j} = 0;
        history_seq.shuffled_songs.(trans_hist_interest{4}){i}{j} = 0;
    
        for ii = 1:length(history_seq.shuffled_songs.syl_order{i}{j})-1
            
            %for transition from ab to ab
            if history_seq.shuffled_songs.syl_order{i}{j}(ii) == tpofinterest{1}(2)...
                    && history_seq.shuffled_songs.syl_order{i}{j}(ii) == history_seq.shuffled_songs.syl_order{i}{j}(ii+1)
                history_seq.shuffled_songs.(trans_hist_interest{1}){i}{j} = ...
                    history_seq.shuffled_songs.(trans_hist_interest{1}){i}{j} + 1;
                history_seq.shuffled_songs.(trans_hist_interest{2}){i}{j} = ...
                    history_seq.shuffled_songs.(trans_hist_interest{2}){i}{j} + 0;
                
                %for transition from ab to ac
            elseif history_seq.shuffled_songs.syl_order{i}{j}(ii) == tpofinterest{1}(2)...
                    && history_seq.shuffled_songs.syl_order{i}{j}(ii) ~= history_seq.shuffled_songs.syl_order{i}{j}(ii+1)
                history_seq.shuffled_songs.(trans_hist_interest{1}){i}{j} = ...
                    history_seq.shuffled_songs.(trans_hist_interest{1}){i}{j} + 0;
                history_seq.shuffled_songs.(trans_hist_interest{2}){i}{j} = ...
                    history_seq.shuffled_songs.(trans_hist_interest{2}){i}{j} + 1;
                
                %for transition from ac to ab
            elseif history_seq.shuffled_songs.syl_order{i}{j}(ii) == tpofinterest{2}(2)...
                    && history_seq.shuffled_songs.syl_order{i}{j}(ii) ~= history_seq.shuffled_songs.syl_order{i}{j}(ii+1)
                history_seq.shuffled_songs.(trans_hist_interest{3}){i}{j} = ...
                    history_seq.shuffled_songs.(trans_hist_interest{1}){i}{j} + 1;
                history_seq.shuffled_songs.(trans_hist_interest{4}){i}{j} = ...
                    history_seq.shuffled_songs.(trans_hist_interest{2}){i}{j} + 0;
                
                %for transition from ac to ac
            elseif history_seq.shuffled_songs.syl_order{i}{j}(ii) == tpofinterest{2}(2)...
                    && history_seq.shuffled_songs.syl_order{i}{j}(ii) == history_seq.shuffled_songs.syl_order{i}{j}(ii+1)
                history_seq.shuffled_songs.(trans_hist_interest{3}){i}{j} = ...
                    history_seq.shuffled_songs.(trans_hist_interest{1}){i}{j} + 0;
                history_seq.shuffled_songs.(trans_hist_interest{4}){i}{j} = ...
                    history_seq.shuffled_songs.(trans_hist_interest{2}){i}{j} + 1;
            end
            
        end
        
    end
    
end

%converts part of transition probabilities from cell to matrix
for i = 1:length(trans_hist_interest)
    for j = 1:length(history_seq.shuffled_songs.(trans_hist_interest{i}))
        history_seq.shuffled_songs.(trans_hist_interest{i}){j} = cell2mat(history_seq.shuffled_songs.(trans_hist_interest{i}){j});
    end
end

%calculates transitional probabilities per song

for i = 1:length(trans_hist_interest)
    for j = 1:length(history_seq.shuffled_songs.(trans_hist_interest{i}))
        if mod(i,2) == 1
            history_seq.shuffled_songs.(trans_hist_interest{i}){j} = history_seq.shuffled_songs.(trans_hist_interest{i}){j}./...
                (history_seq.shuffled_songs.(trans_hist_interest{i}){j}+history_seq.shuffled_songs.(trans_hist_interest{i+1}){j});
        elseif mod(i,2) == 0
            history_seq.shuffled_songs.(trans_hist_interest{i}){j} = ones(size(history_seq.shuffled_songs.(trans_hist_interest{i}){j})) ...
                -history_seq.shuffled_songs.(trans_hist_interest{i-1}){j};
        end
    end    
end

%calculates transitional probabilities per day
for i = 1:length(trans_hist_interest)
    %puts all contents in a giant matrix to calculate transitional
    history_seq.shuffled_songs.(trans_hist_interest{i}) =  cell2mat(history_seq.shuffled_songs.(trans_hist_interest{i}));
    
    %calculates mean per day
    history_seq.shuffled_songs.(trans_hist_interest{i}) = nanmean(history_seq.shuffled_songs.(trans_hist_interest{i}),2);
end


%calculates history depedence (ex: |p(ab|ab)-p(ab|ac)|)

for i = 1:length(tpofinterest)
    switch i
        case 1
            history_seq.shuffled_songs.hist_dep.(tpofinterest{i}).all = abs(history_seq.shuffled_songs.(trans_hist_interest{1})- ...
                history_seq.shuffled_songs.(trans_hist_interest{3}));
        case 2
            history_seq.shuffled_songs.hist_dep.(tpofinterest{i}).all = abs(history_seq.shuffled_songs.(trans_hist_interest{2})- ...
                history_seq.shuffled_songs.(trans_hist_interest{4}));
    end
    
    %calculates 95% CI for shuffled songs
    history_seq.shuffled_songs.hist_dep.(tpofinterest{i}).CI = [prctile(history_seq.shuffled_songs.hist_dep.(tpofinterest{i}).all,5)...
        prctile(history_seq.shuffled_songs.hist_dep.(tpofinterest{i}).all,95)];
end

%% Determines whether bird is history dependent and puts everything in summary structure

if history_seq.hist_dep.(tpofinterest{1}).CI(1) > history_seq.shuffled_songs.hist_dep.(tpofinterest{1}).CI(2)
    history_seq.summary.is_bird_history_dependent = 'yes';
    marker_color = [1 0 0];
elseif history_seq.hist_dep.(tpofinterest{1}).CI(1) <= history_seq.shuffled_songs.hist_dep.(tpofinterest{1}).CI(2)
    history_seq.summary.is_bird_history_dependent = 'no';
    marker_color = [0 0 0];
end

history_seq.summary.median = history_seq.hist_dep.(tpofinterest{1}).median;
history_seq.summary.CI = history_seq.hist_dep.(tpofinterest{1}).CI;
history_seq.summary.shuffled_CI = history_seq.shuffled_songs.hist_dep.(tpofinterest{1}).CI;


%figure for history dependence
figure(1234567), hold on

%line for 95% CI of shuffled songs
line([history_seq.shuffled_songs.hist_dep.(tpofinterest{1}).CI(1) history_seq.shuffled_songs.hist_dep.(tpofinterest{1}).CI(2)],[1 1],...
    'LineWidth', 2, 'Color', [.5 .5 .5])

%line and marker for measured history dependence
plot(history_seq.hist_dep.(tpofinterest{1}).median,1.5,...
    'Marker', 'v', 'MarkerSize', 20, 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', marker_color)
line([history_seq.hist_dep.(tpofinterest{1}).CI(1) history_seq.hist_dep.(tpofinterest{1}).CI(2)], [1.5 1.5],...
    'LineWidth', 2, 'Color', marker_color)
set(gca,'YTick',[])
ylim([0 2])
title(['History dependence of ' bird_name ' on ' date_of_exp])
xlabel(['| p(' tpofinterest{1} '|' tpofinterest{1} ') - p(' tpofinterest{1} '|' tpofinterest{2} ') |'])

%saves figure
saveas(figure(1234567), ['check_sequence_timing_' today_date '/Is_History_Dependent_' bird_name '_' date_of_exp], 'fig')







