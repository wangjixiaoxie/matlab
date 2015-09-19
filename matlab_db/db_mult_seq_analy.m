%To get a quick calculation of transition probabilites with more than 2
%branch points.

%%
dbstop if error


%% Asks for basic parameters
display(cd)
which_computer = input('Bluejay number? ');
multiple_trans.parameters.computer = ['/bluejay' num2str(which_computer) '/dbrady'];
if strcmpi(input(['Is this correct: ' multiple_trans.parameters.computer '?  '],'s'),'y') == 1
else
    multiple_trans.parameters.computer = input('Input path:  ','s');
end
clear which_computer


multiple_trans.parameters.bird_name = input('What is the name of your bird? ', 's');
multiple_trans.parameters.date_of_exp = input('what is the date of your experiment?\n(ex: 03Aug2012)  ', 's');

% A date stamp that will label all your saved variables and figures with the
% time that you did the analysis
multiple_trans.parameters.date_of_analysis = datestr(now, 'ddmmmyyyy_HHMMAM');
multiple_trans.parameters.date_of_analysis = multiple_trans.parameters.date_of_analysis(multiple_trans.parameters.date_of_analysis ~= ' ');

% % Make a directory to store figures and variables
% mkdir(['multi_sequence_' multiple_trans.parameters.date_of_analysis])

%Asks to add a modifying phrase
multiple_trans.parameters.phrase = input('Phrase?  ','s');
if ~isempty(strfind(multiple_trans.parameters.phrase,'''')) == 1
    multiple_trans.parameters.phrase = [];
end

% which batch file do want to analyze (probably batch.catch.keep)
if strcmpi(input('Is your batch file called batch.catch.keep? (y or n) ','s'),'y')
    multiple_trans.parameters.batchfile = 'batch.catch.keep';
else
    multiple_trans.parameters.batchfile = input('what is the name of your batch file?  ', 's');
end

%asks for the number of transition points
multiple_trans.parameters.num_tp = input('How many transition points?  ');

for i = 1:multiple_trans.parameters.num_tp
    %asks for the name of the transition points
    multiple_trans.parameters.tpofinterest{i} = input(['What is the name of transition #' num2str(i) '?  '], 's');
end

if multiple_trans.parameters.num_tp == 2
    if multiple_trans.parameters.tpofinterest{1}(end) == multiple_trans.parameters.tpofinterest{2}(end)
        multiple_trans.parameters.con_or_div = 'con';
    elseif multiple_trans.parameters.tpofinterest{1}(end) ~= multiple_trans.parameters.tpofinterest{2}(end)
        multiple_trans.parameters.con_or_div = 'div';
    end
else
    % asks if it is a convergent or divergent sequence you are interested in
    multiple_trans.parameters.con_or_div = input('Convergence or divergence? (con or div)  ','s');
end

        %Measures the number of syllables from start syllable the syllable of
        %interest is (ex: ab -> 1; abcd -> 3)
for i = 1:multiple_trans.parameters.num_tp
    multiple_trans.parameters.tp_length{i} = length(multiple_trans.parameters.tpofinterest{i})-1;
end


% Make a directory to store figures and variables
mkdir(['multi_sequence_' multiple_trans.parameters.date_of_analysis '_' multiple_trans.parameters.con_or_div '_' multiple_trans.parameters.phrase])
        

%% Gets all syllables per song
multiple_trans.syl_order = getlabels(multiple_trans.parameters.batchfile);

% Gets rid of all syllables that are not the end point of a transition
for i  = 1:length(multiple_trans.syl_order)
    multiple_trans.syl_order{i} = multiple_trans.syl_order{i}(multiple_trans.syl_order{i} ~= '-');
    
    temp_syl_order = [];
    for j = 1:multiple_trans.parameters.num_tp
        if strcmpi(multiple_trans.parameters.con_or_div, 'div')
            temp_syl_order = [temp_syl_order...
                strfind(multiple_trans.syl_order{i}, multiple_trans.parameters.tpofinterest{j})+multiple_trans.parameters.tp_length{j}];
        elseif strcmpi(multiple_trans.parameters.con_or_div, 'con')
            temp_syl_order = [temp_syl_order strfind(multiple_trans.syl_order{i}, multiple_trans.parameters.tpofinterest{j})];
        end
    end
    temp_syl_order = sort(temp_syl_order);
    multiple_trans.syl_order{i} = multiple_trans.syl_order{i}(temp_syl_order);
    clear temp_syl_order

end

% Gets rid of all songs that do not have a transition of interest
multiple_trans.syl_order = multiple_trans.syl_order(~cellfun('isempty',multiple_trans.syl_order));

%% Calculates transition probability

%calculates transition probability each day. This adds the number of each 
%syllable per day and then divides by the
%total number of syllables of interest. The 5th, 50th, and 95th percentiles of
%a permuation test is used to give you the average and CI for the
%transition probability for the day.

%converts syl_order to a matrix. Is a pool of all syllables sung that day
multiple_trans.trans_prob.bootstrap.pool = cell2mat(multiple_trans.syl_order);
multiple_trans.parameters.numbtstrp_trials = 10000;

for j = 1:multiple_trans.parameters.numbtstrp_trials

    %which syllables are being sampled from the pool for this run of the bootstrap procedure
    multiple_trans.trans_prob.bootstrap.sampling{j} =...
        randi(size(multiple_trans.trans_prob.bootstrap.pool,2),[1 size(multiple_trans.trans_prob.bootstrap.pool,2)]);
    
    %the subsequent number of syllables using the sampling above
    multiple_trans.trans_prob.bootstrap.boot_run{j} =...
        multiple_trans.trans_prob.bootstrap.pool(multiple_trans.trans_prob.bootstrap.sampling{j});
    
    %calculates the transition probability for each syllable in the run of
    %the bootstrap procedure
    for i = 1:length(multiple_trans.parameters.tpofinterest)
        %calcualtes the transition probability for syllable of interest for
        %this run of the bootstrap procedure
        if strcmpi(multiple_trans.parameters.con_or_div,'div')
            multiple_trans.trans_prob.bootstrap.(multiple_trans.parameters.tpofinterest{i}){j} =...
                sum(multiple_trans.trans_prob.bootstrap.boot_run{j} == multiple_trans.parameters.tpofinterest{i}(end))...
                ./length(multiple_trans.trans_prob.bootstrap.boot_run{j});
        elseif strcmpi(multiple_trans.parameters.con_or_div,'con')
            multiple_trans.trans_prob.bootstrap.(multiple_trans.parameters.tpofinterest{i}){j} =...
                sum(multiple_trans.trans_prob.bootstrap.boot_run{j} == multiple_trans.parameters.tpofinterest{i}(end-1))...
                ./length(multiple_trans.trans_prob.bootstrap.boot_run{j});
        end
    end

end

%converts the transition probability per day bootstrap data from a cell to
%a matrix (easier to work with)
for i = 1:length(multiple_trans.parameters.tpofinterest)
    multiple_trans.trans_prob.bootstrap.(multiple_trans.parameters.tpofinterest{i}) =...
        cell2mat(multiple_trans.trans_prob.bootstrap.(multiple_trans.parameters.tpofinterest{i}));
end

%gives you the median and CI of the bootstrap procedure
for i = 1:length(multiple_trans.parameters.tpofinterest)
    multiple_trans.trans_prob.median.(multiple_trans.parameters.tpofinterest{i}) =...
        median(multiple_trans.trans_prob.bootstrap.(multiple_trans.parameters.tpofinterest{i}));
    multiple_trans.trans_prob.CI.(multiple_trans.parameters.tpofinterest{i}) = ...
        [prctile(multiple_trans.trans_prob.bootstrap.(multiple_trans.parameters.tpofinterest{i}),2.5)...
        prctile(multiple_trans.trans_prob.bootstrap.(multiple_trans.parameters.tpofinterest{i}),97.5)];
end

%makes histogram of transition probability bootstrap data
for i = 1:length(multiple_trans.parameters.tpofinterest)
    
    %Makes the bins and counts for a histogram (100 bins)
    [multiple_trans.trans_prob.bootstrap.hist.(multiple_trans.parameters.tpofinterest{i}).count...
        multiple_trans.trans_prob.bootstrap.hist.(multiple_trans.parameters.tpofinterest{i}).bins] =...
        hist(multiple_trans.trans_prob.bootstrap.(multiple_trans.parameters.tpofinterest{i}),100);
    
    %Converts the counts into relative frequencies
    multiple_trans.trans_prob.bootstrap.hist.(multiple_trans.parameters.tpofinterest{i}).rel_freq = ...
        multiple_trans.trans_prob.bootstrap.hist.(multiple_trans.parameters.tpofinterest{i}).count./...
        sum(multiple_trans.trans_prob.bootstrap.hist.(multiple_trans.parameters.tpofinterest{i}).count);
    
    
    figure()
    hold on
    bar(multiple_trans.trans_prob.bootstrap.hist.(multiple_trans.parameters.tpofinterest{i}).bins,...
        multiple_trans.trans_prob.bootstrap.hist.(multiple_trans.parameters.tpofinterest{i}).rel_freq)
    db_prepare_bootstrp_hist(...
        ['Transition probability for ' multiple_trans.parameters.tpofinterest{i}...
        ' on ' multiple_trans.parameters.date_of_exp '  (bootstrap' num2str(multiple_trans.parameters.numbtstrp_trials) ' trials)'],...
        'Relative Frequency',...
        'Transition Probability',...
        multiple_trans.trans_prob.CI.(multiple_trans.parameters.tpofinterest{i}),...
        multiple_trans.trans_prob.bootstrap.hist.(multiple_trans.parameters.tpofinterest{i}).rel_freq,...
        multiple_trans.trans_prob.median.(multiple_trans.parameters.tpofinterest{i}))

    saveas(figure(gcf), ['multi_sequence_' multiple_trans.parameters.date_of_analysis '_' multiple_trans.parameters.con_or_div '_'...
        multiple_trans.parameters.phrase ...
        '/Transition_probability_' multiple_trans.parameters.tpofinterest{i} '_' multiple_trans.parameters.date_of_exp], 'fig')
   
end

%% calculates the number of each syllable per trial
multiple_trans.trans_prob.per_trial = zeros(length(multiple_trans.syl_order),length(multiple_trans.parameters.tpofinterest));
for i = 1:length(multiple_trans.syl_order)
    for j = 1:length(multiple_trans.parameters.tpofinterest)
        if strcmpi(multiple_trans.parameters.con_or_div,'div')
            multiple_trans.trans_prob.per_trial(i,j) = sum(multiple_trans.syl_order{i} == multiple_trans.parameters.tpofinterest{j}(end));
        elseif strcmpi(multiple_trans.parameters.con_or_div,'con')
            multiple_trans.trans_prob.per_trial(i,j) = sum(multiple_trans.syl_order{i} == multiple_trans.parameters.tpofinterest{j}(end-1));
        end
    end
end

%% calculates the relative frequency of each syllable transistion per song
for i = 1:size(multiple_trans.trans_prob.per_trial,2)
    for j = 1:size(multiple_trans.trans_prob.per_trial,1)
        multiple_trans.trans_prob.rel_freq.rel_freq(j,i) = multiple_trans.trans_prob.per_trial(j,i)./sum(multiple_trans.trans_prob.per_trial(j,:));
    end
end

%does a little bootstrapping to find and median, std, and cv of relative frequency
%throughout the day

multiple_trans.trans_prob.rel_freq.bootstrap =...
    bootstrp(multiple_trans.parameters.numbtstrp_trials,@(y) [median(y) std(y)],multiple_trans.trans_prob.rel_freq.rel_freq);
multiple_trans.trans_prob.rel_freq.boot_CI =...
    bootci(multiple_trans.parameters.numbtstrp_trials, @(y) [median(y) std(y)], multiple_trans.trans_prob.rel_freq.rel_freq);

for i = 1:length(multiple_trans.parameters.tpofinterest)
    %median relative freq
    multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{i}).median =...
        median(multiple_trans.trans_prob.rel_freq.bootstrap(:,i));
    multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{i}).median_CI = ...
        multiple_trans.trans_prob.rel_freq.boot_CI(:,i);
    
    %std relative freq
    multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{i}).sd = ...
        median(multiple_trans.trans_prob.rel_freq.bootstrap(:,length(multiple_trans.parameters.tpofinterest)+i));
    multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{i}).sd_CI = ...
        multiple_trans.trans_prob.rel_freq.boot_CI(:,length(multiple_trans.parameters.tpofinterest)+i);
    
%     %cv of relative freq
%     multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{i}).cv = ...
%         nanmedian(multiple_trans.trans_prob.rel_freq.bootstrap(:,2*length(multiple_trans.parameters.tpofinterest)+i));
%     multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{i}).cv_CI = ...
%         multiple_trans.trans_prob.rel_freq.boot_CI(:,2*length(multiple_trans.parameters.tpofinterest)+i);
end

%% To see if there is a change in transition probability over the day: Part 1 (every song)

%calcualtes the transition probability per song throughout the day
for j = 1:size(multiple_trans.trans_prob.per_trial,1)
    for i = 1:length(multiple_trans.parameters.tpofinterest) 
        multiple_trans.trans_prob.day.trans_prob_per_song(j,i) =...
            multiple_trans.trans_prob.per_trial(j,i)./sum(multiple_trans.trans_prob.per_trial(j,:),2);    
    end
end

%makes figure of transition probability per song throughout the day
figure(), hold on
for i = 1:length(multiple_trans.parameters.tpofinterest)
    plot(multiple_trans.trans_prob.day.trans_prob_per_song(:,i),...
        '-o',...
        'MarkerSize', 10,...
        'Color', [1-((length(multiple_trans.parameters.tpofinterest)-i)./length(multiple_trans.parameters.tpofinterest))...
        0 0+((length(multiple_trans.parameters.tpofinterest)-i)./length(multiple_trans.parameters.tpofinterest))])
end
legend(multiple_trans.parameters.tpofinterest)
title(['Tranisition probability of each song throughout ' multiple_trans.parameters.date_of_exp])
xlabel('Song Number')
ylabel('Transition probability')
ylim([-0.1 1.1])
saveas(figure(gcf), ['multi_sequence_' multiple_trans.parameters.date_of_analysis '_' multiple_trans.parameters.con_or_div '_'...
    multiple_trans.parameters.phrase ...
    '/Trans_prob_per_song_' multiple_trans.parameters.date_of_exp], 'fig')

%% To see if there is a change in transition probability over the day: Part 2 (cumulative average, across syllables)

%Calculates cumulative sum of a certain syllable over the pool of syllables
for j = 1:length(multiple_trans.parameters.tpofinterest)
    if strcmpi(multiple_trans.parameters.con_or_div,'div')
        multiple_trans.trans_prob.day.cumulative_trans_prob.(multiple_trans.parameters.tpofinterest{j}) =...
            cumsum(multiple_trans.trans_prob.bootstrap.pool == multiple_trans.parameters.tpofinterest{j}(end));
    elseif strcmpi(multiple_trans.parameters.con_or_div,'con')
        multiple_trans.trans_prob.day.cumulative_trans_prob.(multiple_trans.parameters.tpofinterest{j}) =...
            cumsum(multiple_trans.trans_prob.bootstrap.pool == multiple_trans.parameters.tpofinterest{j}(end-1));
    end
end

%Calculates cumulative fraction of a certain syllable over the day
for j = 1:length(multiple_trans.parameters.tpofinterest)
    for i = 1:length(multiple_trans.trans_prob.day.cumulative_trans_prob.(multiple_trans.parameters.tpofinterest{j}))
        multiple_trans.trans_prob.day.cumulative_trans_prob.(multiple_trans.parameters.tpofinterest{j})(i) = ...
            multiple_trans.trans_prob.day.cumulative_trans_prob.(multiple_trans.parameters.tpofinterest{j})(i)./i;
    end
end

%makes figure of cumulative transition probability
figure(), hold on
for i = 1:length(multiple_trans.parameters.tpofinterest)
    plot(multiple_trans.trans_prob.day.cumulative_trans_prob.(multiple_trans.parameters.tpofinterest{i}),...
        '-o',...
        'MarkerSize', 6,...
        'Color', [1-((length(multiple_trans.parameters.tpofinterest)-i)./length(multiple_trans.parameters.tpofinterest))...
        0 0+((length(multiple_trans.parameters.tpofinterest)-i)./length(multiple_trans.parameters.tpofinterest))])
end
legend(multiple_trans.parameters.tpofinterest)
title(['Cumulative tranisition probability on ' multiple_trans.parameters.date_of_exp])
xlabel('Syllable Number')
ylabel('Transition probability')
ylim([-0.1 1.1])
saveas(figure(gcf), ['multi_sequence_' multiple_trans.parameters.date_of_analysis '_' multiple_trans.parameters.con_or_div '_'...
    multiple_trans.parameters.phrase ...
    '/Cumul_trans_prob_' multiple_trans.parameters.date_of_exp], 'fig')

%% To see if there is a change in transition probability over the day: Part 3 (running average, across songs)

multiple_trans.parameters.running_average_window = 5;

for j = 1:length(multiple_trans.parameters.tpofinterest)
    multiple_trans.trans_prob.day.run_avg_trans_prob(:,j) =...
        db_runningaverage(multiple_trans.trans_prob.day.trans_prob_per_song(:,j),...
        multiple_trans.parameters.running_average_window,1);
end

%makes figure of running average transition probability
figure(), hold on
for i = 1:length(multiple_trans.parameters.tpofinterest)
    plot(multiple_trans.trans_prob.day.run_avg_trans_prob(:,i),...
        '-o',...
        'MarkerSize', 10,...
        'Color', [1-((length(multiple_trans.parameters.tpofinterest)-i)./length(multiple_trans.parameters.tpofinterest))...
        0 0+((length(multiple_trans.parameters.tpofinterest)-i)./length(multiple_trans.parameters.tpofinterest))])
end
legend(multiple_trans.parameters.tpofinterest)
title(['Running average tranisition probability on ' multiple_trans.parameters.date_of_exp...
    '   (window size: ' num2str(multiple_trans.parameters.running_average_window) ')'])
xlabel('Song Number')
ylabel('Transition probability')
ylim([-0.1 1.1])
saveas(figure(gcf), ['multi_sequence_' multiple_trans.parameters.date_of_analysis '_' multiple_trans.parameters.con_or_div '_'...
    multiple_trans.parameters.phrase ...
    '/Run_avg_trans_prob_' multiple_trans.parameters.date_of_exp], 'fig')

%% Calculates entropy ( Sigma -p log2 p)

for j = 1:length(multiple_trans.parameters.tpofinterest)
    multiple_trans.entropy.(multiple_trans.parameters.tpofinterest{j}) =...
        -1.*multiple_trans.trans_prob.bootstrap.(multiple_trans.parameters.tpofinterest{j}).*log2(...
        multiple_trans.trans_prob.bootstrap.(multiple_trans.parameters.tpofinterest{j}));
    %converts -p log2 p to 0 from NaN if p is 0 (the limit)
    multiple_trans.entropy.(multiple_trans.parameters.tpofinterest{j})(isnan(multiple_trans.entropy.(multiple_trans.parameters.tpofinterest{j})))...
        = 0;
end

multiple_trans.entropy.entropy = [];
for j = 1:length(multiple_trans.parameters.tpofinterest)
    multiple_trans.entropy.entropy(j,:) = multiple_trans.entropy.(multiple_trans.parameters.tpofinterest{j});
end
multiple_trans.entropy.entropy = sum(multiple_trans.entropy.entropy);


multiple_trans.entropy.median = nanmedian(multiple_trans.entropy.entropy);
multiple_trans.entropy.CI = [prctile(multiple_trans.entropy.entropy,2.5) prctile(multiple_trans.entropy.entropy,97.5)];

%Makes histogram of entropy
%Makes the bins and counts for a histogram (100 bins)
[multiple_trans.entropy.hist.count multiple_trans.entropy.hist.bins] = hist(multiple_trans.entropy.entropy,100);

%Converts the counts into relative frequencies
multiple_trans.entropy.hist.rel_freq = multiple_trans.entropy.hist.count ./ sum(multiple_trans.entropy.hist.count);


figure()
hold on
bar(multiple_trans.entropy.hist.bins, multiple_trans.entropy.hist.rel_freq)
db_prepare_bootstrp_hist(...
    ['Entropy for ' multiple_trans.parameters.bird_name ' on '...
    multiple_trans.parameters.date_of_exp '  (bootstrap' num2str(multiple_trans.parameters.numbtstrp_trials) ' trials)'],...
    'Relative Frequency',...
    'Entropy (Sigma -p log2 p)',...
    multiple_trans.entropy.CI,...
    multiple_trans.entropy.hist.rel_freq,...
    multiple_trans.entropy.median)

saveas(figure(gcf), ['multi_sequence_' multiple_trans.parameters.date_of_analysis '_' multiple_trans.parameters.con_or_div '_' ...
    multiple_trans.parameters.phrase '/Entropy_'...
    multiple_trans.parameters.bird_name '_' multiple_trans.parameters.date_of_exp], 'fig')



%% Calculating History Dependence with real data (| p(ab|ab) - p(ab|ac) - p(ab|ad) |)

for i = 1:length(multiple_trans.parameters.tpofinterest)
    
    %Makes a cell with other possible transitions besides the one of
    %current interest
    multiple_trans.hist_dep.transitions =...
        multiple_trans.parameters.tpofinterest;
    
    %Calculates the number of a certain transition per song (ex: ab|ab or
    %ac|ad)
    for j = 1:length(multiple_trans.hist_dep.transitions)
        for k = 1:length(multiple_trans.syl_order)
            if strcmpi(multiple_trans.parameters.con_or_div,'div')
                multiple_trans.hist_dep.count.([multiple_trans.parameters.tpofinterest{i} '_to_' multiple_trans.hist_dep.transitions{j}])(k) =...
                    length(strfind(multiple_trans.syl_order{k}, [multiple_trans.parameters.tpofinterest{i}(end) ...
                    multiple_trans.hist_dep.transitions{j}(end)]));
            elseif strcmpi(multiple_trans.parameters.con_or_div,'con')
                multiple_trans.hist_dep.count.([multiple_trans.parameters.tpofinterest{i} '_to_' multiple_trans.hist_dep.transitions{j}])(k) =...
                    length(strfind(multiple_trans.syl_order{k}, [multiple_trans.parameters.tpofinterest{i}(end-1) ...
                    multiple_trans.hist_dep.transitions{j}(end-1)]));
            end
        end
    end
    
end

%Makes a matrix with the names of all the possible history dep comparisons
multiple_trans.hist_dep.possible_history_seq = fieldnames(multiple_trans.hist_dep.count);
multiple_trans.hist_dep.possible_history_seq = cell2mat(multiple_trans.hist_dep.possible_history_seq);

%Calculates the relative frequency of a history dependent transition per
%trial
for i = 1:size(multiple_trans.hist_dep.possible_history_seq,1)
    for j = 1:length(multiple_trans.syl_order)
        
        total_this_song = [];
        if strcmpi(multiple_trans.parameters.con_or_div, 'div')
            which_comparisons = find(multiple_trans.hist_dep.possible_history_seq(:,length(multiple_trans.parameters.tpofinterest{1}))...
                == multiple_trans.hist_dep.possible_history_seq(i,length(multiple_trans.parameters.tpofinterest{1})));
        elseif strcmpi(multiple_trans.parameters.con_or_div, 'con')
            which_comparisons = find(multiple_trans.hist_dep.possible_history_seq(:,length(multiple_trans.parameters.tpofinterest{1})-1)...
            == multiple_trans.hist_dep.possible_history_seq(i,length(multiple_trans.parameters.tpofinterest{1})-1));
        end
        for k = 1:length(which_comparisons)
            total_this_song = [total_this_song multiple_trans.hist_dep.count.(multiple_trans.hist_dep.possible_history_seq(which_comparisons(k),:))(j)];
        end
        total_this_song = sum(total_this_song);
        which_comparisons = [];
        
        multiple_trans.hist_dep.rel_freq.(multiple_trans.hist_dep.possible_history_seq(i,:))(j) = ...
            multiple_trans.hist_dep.count.(multiple_trans.hist_dep.possible_history_seq(i,:))(j)...
            ./total_this_song;
    end
end

%deletes temporary variables
clear which_comparisons
clear total_this_song

%Gets rid of NaN
for i = 1:size(multiple_trans.hist_dep.possible_history_seq,1)
    multiple_trans.hist_dep.rel_freq.(multiple_trans.hist_dep.possible_history_seq(i,:)) = ...
        multiple_trans.hist_dep.rel_freq.(multiple_trans.hist_dep.possible_history_seq(i,:))...
        (~isnan(multiple_trans.hist_dep.rel_freq.(multiple_trans.hist_dep.possible_history_seq(i,:))));
    
    if isempty(multiple_trans.hist_dep.rel_freq.(multiple_trans.hist_dep.possible_history_seq(i,:))) == 1
        multiple_trans.hist_dep.rel_freq.(multiple_trans.hist_dep.possible_history_seq(i,:)) = 0;
    else
    end
end

%%%Calculates history dependence using a bootstrap procedure
%First does bootstrap for certain history dependent transition
for i = 1:size(multiple_trans.hist_dep.possible_history_seq,1)
    for k = 1:multiple_trans.parameters.numbtstrp_trials
        multiple_trans.hist_dep.bootstrp.(multiple_trans.hist_dep.possible_history_seq(i,:)).sampling{k} = ...
            randi(size(multiple_trans.hist_dep.rel_freq.(multiple_trans.hist_dep.possible_history_seq(i,:)),2),...
            [1 size(multiple_trans.hist_dep.rel_freq.(multiple_trans.hist_dep.possible_history_seq(i,:)),2)]);
        
        multiple_trans.hist_dep.bootstrp.(multiple_trans.hist_dep.possible_history_seq(i,:)).boot_run{k} = ...
            multiple_trans.hist_dep.rel_freq.(multiple_trans.hist_dep.possible_history_seq(i,:))...
            (multiple_trans.hist_dep.bootstrp.(multiple_trans.hist_dep.possible_history_seq(i,:)).sampling{k});
        
        multiple_trans.hist_dep.bootstrp.(multiple_trans.hist_dep.possible_history_seq(i,:)).mean(k) = ...
            mean(multiple_trans.hist_dep.bootstrp.(multiple_trans.hist_dep.possible_history_seq(i,:)).boot_run{k});
    end
end

%Calculates history dependence
for i = 1:length(multiple_trans.parameters.tpofinterest)
    
    if strcmpi(multiple_trans.parameters.con_or_div,'div')
        which_comparisons = find(multiple_trans.hist_dep.possible_history_seq(:,end)...
            == multiple_trans.hist_dep.possible_history_seq(i,end));
    elseif strcmpi(multiple_trans.parameters.con_or_div,'con')
        which_comparisons = find(multiple_trans.hist_dep.possible_history_seq(:,end-1)...
            == multiple_trans.hist_dep.possible_history_seq(i,end-1));
    end
    
    diff_matrix = [];
    
    for j = 1:length(which_comparisons)
        diff_matrix = [diff_matrix;...
            multiple_trans.hist_dep.bootstrp.(multiple_trans.hist_dep.possible_history_seq(which_comparisons(j),:)).mean];
    end
    
    diff_matrix = flipdim(sortrows(diff_matrix),1);
    
    for k = 1:size(diff_matrix,1)-1
        diff_matrix(k+1,:) = diff_matrix(k,:)-diff_matrix(k+1,:);
        diff_matrix(k,:) = 0;
    end
    
    diff_matrix = sum(abs(diff_matrix));
    
    multiple_trans.hist_dep.bootstrp.(multiple_trans.parameters.tpofinterest{i}) = diff_matrix;
end

clear diff_matrix
clear which_comparisons

%Calculates the 5th, 50th, and 95th percentiles of history dependence
for i = 1:length(multiple_trans.parameters.tpofinterest)
    %median of bootstrap distro
    multiple_trans.hist_dep.median.(multiple_trans.parameters.tpofinterest{i}) =...
        median(multiple_trans.hist_dep.bootstrp.(multiple_trans.parameters.tpofinterest{i}));
    
    %CI of bootstrap distro
    multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{i}) = ...
        [prctile(multiple_trans.hist_dep.bootstrp.(multiple_trans.parameters.tpofinterest{i}),2.5)...
        prctile(multiple_trans.hist_dep.bootstrp.(multiple_trans.parameters.tpofinterest{i}),97.5)];
end

%makes histogram of history dependence bootstrap data
for i = 1:length(multiple_trans.parameters.tpofinterest)
    
    %Makes the bins and counts for a histogram (100 bins)
    [multiple_trans.hist_dep.hist.(multiple_trans.parameters.tpofinterest{i}).count...
        multiple_trans.hist_dep.hist.(multiple_trans.parameters.tpofinterest{i}).bins] =...
        hist(multiple_trans.hist_dep.bootstrp.(multiple_trans.parameters.tpofinterest{i}),100);
    
    %Converts the counts into relative frequencies
    multiple_trans.hist_dep.hist.(multiple_trans.parameters.tpofinterest{i}).rel_freq = ...
        multiple_trans.hist_dep.hist.(multiple_trans.parameters.tpofinterest{i}).count./...
        sum(multiple_trans.hist_dep.hist.(multiple_trans.parameters.tpofinterest{i}).count);
    
    
    figure()
    hold on
    bar(multiple_trans.hist_dep.hist.(multiple_trans.parameters.tpofinterest{i}).bins,...
        multiple_trans.hist_dep.hist.(multiple_trans.parameters.tpofinterest{i}).rel_freq)
    db_prepare_bootstrp_hist(...
        ['History dependence for ' multiple_trans.parameters.tpofinterest{i}...
        ' on ' multiple_trans.parameters.date_of_exp '  (bootstrap' num2str(multiple_trans.parameters.numbtstrp_trials) ' trials)'],...
        'Relative Frequency',...
        'History Dependence',...
        multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{i}),...
        multiple_trans.hist_dep.hist.(multiple_trans.parameters.tpofinterest{i}).rel_freq,...
        multiple_trans.hist_dep.median.(multiple_trans.parameters.tpofinterest{i}))

    saveas(figure(gcf), ['multi_sequence_' multiple_trans.parameters.date_of_analysis '_' multiple_trans.parameters.con_or_div '_'...
        multiple_trans.parameters.phrase '/History_dependence_'...
        multiple_trans.parameters.tpofinterest{i} '_' multiple_trans.parameters.date_of_exp], 'fig')
   
end

%% Calculating shuffled history dependence - across songs shuffle

% %Combines all syllables sung in a day into a pool and shuffles them
% for i = 1:multiple_trans.parameters.numbtstrp_trials
%     multiple_trans.shuffled_songs.pool(i,:) = multiple_trans.trans_prob.bootstrap.pool(randperm(length(multiple_trans.trans_prob.bootstrap.pool)));
% end
% 
% %Finds the lengths of the individual songs sung throughout the day
% multiple_trans.shuffled_songs.lengths_of_songs = 0;
% for i = 1:length(multiple_trans.syl_order)
%     multiple_trans.shuffled_songs.lengths_of_songs(i+1) = length(multiple_trans.syl_order{i});
% end
% multiple_trans.shuffled_songs.lengths_of_songs = cumsum(multiple_trans.shuffled_songs.lengths_of_songs);
% 
% %Reconstitutes the songs
% for i = 1:length(multiple_trans.syl_order)
%     multiple_trans.shuffled_songs.songs{i} = multiple_trans.shuffled_songs.pool(:,multiple_trans.shuffled_songs.lengths_of_songs(i)+1:...
%         multiple_trans.shuffled_songs.lengths_of_songs(i+1));
% end
% 
% 
% %Makes a cell with other possible transitions besides the one of
% %current interest
% multiple_trans.shuffled_songs.transitions =...
%     multiple_trans.parameters.tpofinterest;
% 
% 
% for i = 1:length(multiple_trans.parameters.tpofinterest)
%     
%     %Calculates the number of a certain transition per song (ex: ab|ab or
%     %ac|ad)
%     for j = 1:length(multiple_trans.hist_dep.transitions)
%         for k = 1:length(multiple_trans.shuffled_songs.songs)
%             for m = 1:multiple_trans.parameters.numbtstrp_trials
%                 multiple_trans.shuffled_songs.count.([multiple_trans.parameters.tpofinterest{i} '_to_' multiple_trans.hist_dep.transitions{j}]){k}(m) =...
%                     length(strfind(multiple_trans.shuffled_songs.songs{k}(m,:), [multiple_trans.parameters.tpofinterest{i}(end) ...
%                     multiple_trans.hist_dep.transitions{j}(end)]));
%             end
%         end
%     end
%     
% end
% 
% % Makes a matrix with the names of all the possible history dep comparisons
% multiple_trans.shuffled_songs.possible_history_seq = fieldnames(multiple_trans.shuffled_songs.count);
% multiple_trans.shuffled_songs.possible_history_seq = cell2mat(multiple_trans.shuffled_songs.possible_history_seq);
% 
% %Calculates the relative frequency of a history dependent transition per
% %trial
% for i = 1:size(multiple_trans.shuffled_songs.possible_history_seq,1)
%     for j = 1:length(multiple_trans.shuffled_songs.songs)
%         
%             total_this_song = [];
%             which_comparisons = find(multiple_trans.shuffled_songs.possible_history_seq(:,length(multiple_trans.parameters.tpofinterest{1}))...
%                 == multiple_trans.shuffled_songs.possible_history_seq(i,length(multiple_trans.parameters.tpofinterest{1})));
%             for k = 1:length(which_comparisons)
%                 total_this_song = [total_this_song;...
%                     multiple_trans.shuffled_songs.count.(multiple_trans.shuffled_songs.possible_history_seq(which_comparisons(k),:)){j}];
%             end
%             total_this_song = sum(total_this_song,1);
%             which_comparisons = [];
%             
%             multiple_trans.shuffled_songs.rel_freq.(multiple_trans.shuffled_songs.possible_history_seq(i,:)){j} = ...
%                 multiple_trans.shuffled_songs.count.(multiple_trans.shuffled_songs.possible_history_seq(i,:)){j}...
%                 ./total_this_song;
%             
% %             multiple_trans.shuffled_songs.rel_freq.(multiple_trans.shuffled_songs.possible_history_seq(i,:)){j} = ...
% %                 multiple_trans.shuffled_songs.rel_freq.(multiple_trans.shuffled_songs.possible_history_seq(i,:)){j}(~isnan(...
% %                 multiple_trans.shuffled_songs.rel_freq.(multiple_trans.shuffled_songs.possible_history_seq(i,:)){j}));
%     end    
%     
%     %organizes the trials into a pseudo day list
%     multiple_trans.shuffled_songs.rel_freq.(multiple_trans.shuffled_songs.possible_history_seq(i,:)) =...
%         cell2mat(multiple_trans.shuffled_songs.rel_freq.(multiple_trans.shuffled_songs.possible_history_seq(i,:)));
%     multiple_trans.shuffled_songs.rel_freq.(multiple_trans.shuffled_songs.possible_history_seq(i,:)) = reshape(...
%         multiple_trans.shuffled_songs.rel_freq.(multiple_trans.shuffled_songs.possible_history_seq(i,:)),...
%         length(multiple_trans.shuffled_songs.songs),multiple_trans.parameters.numbtstrp_trials);
%     
%     %calculates the nanmean for each pseudo day and gets rid of NaN trials
%     multiple_trans.shuffled_songs.mean_each_day.(multiple_trans.shuffled_songs.possible_history_seq(i,:)) = ...
%         nanmean(multiple_trans.shuffled_songs.rel_freq.(multiple_trans.shuffled_songs.possible_history_seq(i,:)),1);
%     multiple_trans.shuffled_songs.mean_each_day.(multiple_trans.shuffled_songs.possible_history_seq(i,:)) = ...
%         multiple_trans.shuffled_songs.mean_each_day.(multiple_trans.shuffled_songs.possible_history_seq(i,:))(~isnan(...
%         multiple_trans.shuffled_songs.mean_each_day.(multiple_trans.shuffled_songs.possible_history_seq(i,:))));
%     
% end
% 
% %deletes temporary variables
% clear which_comparisons
% clear total_this_song
% 
% 
% %Calculates history dependence
% for i = 1:length(multiple_trans.parameters.tpofinterest)
%     
%     which_comparisons = find(multiple_trans.shuffled_songs.possible_history_seq(:,end)...
%         == multiple_trans.shuffled_songs.possible_history_seq(i,end));
%     diff_matrix = [];
%     
%     for j = 1:length(which_comparisons)
%         diff_matrix = [diff_matrix;...
%             multiple_trans.shuffled_songs.mean_each_day.(multiple_trans.shuffled_songs.possible_history_seq(which_comparisons(j),:))];
%     end
%     
%     diff_matrix = abs(diff(diff_matrix));
%     multiple_trans.shuffled_songs.(multiple_trans.parameters.tpofinterest{i}) = diff_matrix;
% end
% 
% clear diff_matrix
% clear which_comparisons
% 
% %Calculates 95% CI for shuffled history depedence
% 
% for i = 1:length(multiple_trans.parameters.tpofinterest)
%     
%     multiple_trans.shuffled_songs.CI.(multiple_trans.parameters.tpofinterest{i}) = [...
%         prctile(multiple_trans.shuffled_songs.(multiple_trans.parameters.tpofinterest{i}),2.5)... 
%         prctile(multiple_trans.shuffled_songs.(multiple_trans.parameters.tpofinterest{i}),97.5)];
% 
% end
% 
% %% Graph showing CI of history dependence and shuffled history dependence
% 
% 
% for i = 1:length(multiple_trans.parameters.tpofinterest)
%     if multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{i})(1) >...
%             multiple_trans.shuffled_songs.CI.(multiple_trans.parameters.tpofinterest{i})(2)
%         multiple_trans.summary.is_bird_history_dependent = 'yes';
%         marker_color = [1 0 0];
%     elseif multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{i})(1) <=...
%             multiple_trans.shuffled_songs.CI.(multiple_trans.parameters.tpofinterest{i})(2)
%         multiple_trans.summary.is_bird_history_dependent = 'no';
%         marker_color = [0 0 0];
%     end
%     
%     
%     %figure for history dependence
%     figure(1234567+i), hold on
%     
%     %line for 95% CI of shuffled songs
%     line([multiple_trans.shuffled_songs.CI.(multiple_trans.parameters.tpofinterest{i})(1)...
%         multiple_trans.shuffled_songs.CI.(multiple_trans.parameters.tpofinterest{i})(2)],...
%         [1 1],...
%         'LineWidth', 2, 'Color', [.5 .5 .5])
%     
%     %line and marker for measured history dependence
%     plot( multiple_trans.hist_dep.median.(multiple_trans.parameters.tpofinterest{i}),1.5,...
%         'Marker', 'v', 'MarkerSize', 20, 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', marker_color)
%     line([multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{i})(1)...
%         multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{i})(2)], [1.5 1.5],...
%         'LineWidth', 2, 'Color', marker_color)
%     set(gca,'YTick',[])
%     ylim([0 2])
%     title(['History dependence of ' multiple_trans.parameters.bird_name...
%         ' for ' multiple_trans.parameters.tpofinterest{i} ' on ' multiple_trans.parameters.date_of_exp...
%         '    (across song shuffle)'])
%     xlabel(['History dependence (ex: | p(ab|ab) - p(ab|ac) - p(ab|ad) |)'])
%     
%     %saves figure
%     saveas(figure(1234567+i), ['multi_sequence_' multiple_trans.parameters.date_of_analysis...
%         '/Is_History_Dependent_(across_shuffle)_' multiple_trans.parameters.tpofinterest{i} '_' multiple_trans.parameters.date_of_exp], 'fig')
% end



%% %% Calculate shuffled history depedence - within song shuffle

%shuffles each song numbtstrp_trials times
for i = 1:length(multiple_trans.syl_order)
    if length(multiple_trans.syl_order{i}) > 1
        for j = 1:multiple_trans.parameters.numbtstrp_trials
            multiple_trans.within_shuffled.songs{i}(j,:) = multiple_trans.syl_order{i}...
                (randperm(length(multiple_trans.syl_order{i})));
        end
    else
    end
end

%gets rid of empty cells (songs with no transitions)
multiple_trans.within_shuffled.songs = multiple_trans.within_shuffled.songs(~cellfun('isempty',multiple_trans.within_shuffled.songs));

%Makes a cell with other possible transitions besides the one of
%current interest
multiple_trans.within_shuffled.transitions =...
    multiple_trans.parameters.tpofinterest;

%Calculates the number of a certain transition per song (ex: ab|ab or ac|ad)
for i = 1:length(multiple_trans.parameters.tpofinterest)
    

    for j = 1:length(multiple_trans.within_shuffled.transitions)
        for k = 1:length(multiple_trans.within_shuffled.songs)
            for m = 1:multiple_trans.parameters.numbtstrp_trials
                if strcmpi(multiple_trans.parameters.con_or_div,'div')
                    multiple_trans.within_shuffled.count.([multiple_trans.parameters.tpofinterest{i} '_to_' multiple_trans.hist_dep.transitions{j}]){k}(m) =...
                        length(strfind(multiple_trans.within_shuffled.songs{k}(m,:), [multiple_trans.parameters.tpofinterest{i}(end) ...
                        multiple_trans.hist_dep.transitions{j}(end)]));
                elseif strcmpi(multiple_trans.parameters.con_or_div,'con')
                    multiple_trans.within_shuffled.count.([multiple_trans.parameters.tpofinterest{i} '_to_' multiple_trans.hist_dep.transitions{j}]){k}(m) =...
                        length(strfind(multiple_trans.within_shuffled.songs{k}(m,:), [multiple_trans.parameters.tpofinterest{i}(end-1) ...
                        multiple_trans.hist_dep.transitions{j}(end-1)]));
                end
            end
        end
    end
    
end


%Makes a matrix with the names of all the possible history dep comparisons
multiple_trans.within_shuffled.possible_history_seq = fieldnames(multiple_trans.within_shuffled.count);
multiple_trans.within_shuffled.possible_history_seq = cell2mat(multiple_trans.within_shuffled.possible_history_seq);

%Calculates the relative frequency of a history dependent transition per
%trial
for i = 1:size(multiple_trans.within_shuffled.possible_history_seq,1)
    for j = 1:length(multiple_trans.within_shuffled.songs)
        
            total_this_song = [];
            if strcmpi(multiple_trans.parameters.con_or_div, 'div')
                which_comparisons = find(multiple_trans.within_shuffled.possible_history_seq(:,length(multiple_trans.parameters.tpofinterest{1}))...
                    == multiple_trans.within_shuffled.possible_history_seq(i,length(multiple_trans.parameters.tpofinterest{1})));
            elseif strcmpi(multiple_trans.parameters.con_or_div, 'con')
                which_comparisons = find(multiple_trans.within_shuffled.possible_history_seq(:,length(multiple_trans.parameters.tpofinterest{1})-1)...
                    == multiple_trans.within_shuffled.possible_history_seq(i,length(multiple_trans.parameters.tpofinterest{1})-1));
            end

            for k = 1:length(which_comparisons)
                total_this_song = [total_this_song;...
                    multiple_trans.within_shuffled.count.(multiple_trans.within_shuffled.possible_history_seq(which_comparisons(k),:)){j}];
            end
            total_this_song = sum(total_this_song,1);
            which_comparisons = [];
            
            multiple_trans.within_shuffled.rel_freq.(multiple_trans.within_shuffled.possible_history_seq(i,:))(j,:) = ...
                multiple_trans.within_shuffled.count.(multiple_trans.within_shuffled.possible_history_seq(i,:)){j}...
                ./total_this_song;
            
%             for ii = 1:length(multiple_trans.within_shuffled.rel_freq.(multiple_trans.within_shuffled.possible_history_seq(i,:))(j,:))
%                 if isnan(multiple_trans.within_shuffled.rel_freq.(multiple_trans.within_shuffled.possible_history_seq(i,:))(j,ii)) == 1
%                     multiple_trans.within_shuffled.rel_freq.(multiple_trans.within_shuffled.possible_history_seq(i,:))(j,ii) = 0;
%                 else
%                 end
%             end
    end    
    
%     %organizes the trials into a pseudo day list
%     multiple_trans.within_shuffled.rel_freq.(multiple_trans.within_shuffled.possible_history_seq(i,:)) =...
%         cell2mat(multiple_trans.within_shuffled.rel_freq.(multiple_trans.within_shuffled.possible_history_seq(i,:)));
%     multiple_trans.within_shuffled.rel_freq.(multiple_trans.within_shuffled.possible_history_seq(i,:)) = reshape(...
%         multiple_trans.within_shuffled.rel_freq.(multiple_trans.within_shuffled.possible_history_seq(i,:)),...
%         length(multiple_trans.within_shuffled.songs),multiple_trans.parameters.numbtstrp_trials);
    
    %calculates the nanmean for each pseudo day and gets rid of NaN trials
    multiple_trans.within_shuffled.mean_each_day.(multiple_trans.within_shuffled.possible_history_seq(i,:)) = ...
        nanmean(multiple_trans.within_shuffled.rel_freq.(multiple_trans.within_shuffled.possible_history_seq(i,:)),1);
%     multiple_trans.within_shuffled.mean_each_day.(multiple_trans.within_shuffled.possible_history_seq(i,:)) = ...
%         multiple_trans.within_shuffled.mean_each_day.(multiple_trans.within_shuffled.possible_history_seq(i,:))(~isnan(...
%         multiple_trans.within_shuffled.mean_each_day.(multiple_trans.within_shuffled.possible_history_seq(i,:))));
    
end

%deletes temporary variables
clear which_comparisons
clear total_this_song


%Calculates history dependence
for i = 1:length(multiple_trans.parameters.tpofinterest)
    
    if strcmpi(multiple_trans.parameters.con_or_div,'div')
        which_comparisons = find(multiple_trans.within_shuffled.possible_history_seq(:,end)...
            == multiple_trans.within_shuffled.possible_history_seq(i,end));
    elseif strcmpi(multiple_trans.parameters.con_or_div,'con')
        which_comparisons = find(multiple_trans.within_shuffled.possible_history_seq(:,end-1)...
            == multiple_trans.within_shuffled.possible_history_seq(i,end-1));
    end
    diff_matrix = [];
    
    for j = 1:length(which_comparisons)
        diff_matrix = [diff_matrix;...
            multiple_trans.within_shuffled.mean_each_day.(multiple_trans.within_shuffled.possible_history_seq(which_comparisons(j),:))];
    end
    
    diff_matrix = flipdim(sortrows(diff_matrix),1);
    
    for k = 1:size(diff_matrix,1)-1
        diff_matrix(k+1,:) = diff_matrix(k,:)-diff_matrix(k+1,:);
        diff_matrix(k,:) = 0;
    end
    
    diff_matrix = sum(abs(diff_matrix));
    
    multiple_trans.within_shuffled.(multiple_trans.parameters.tpofinterest{i}) = diff_matrix;
end

clear diff_matrix
clear which_comparisons

%Calculates 95% CI for shuffled history depedence

for i = 1:length(multiple_trans.parameters.tpofinterest)
    
    multiple_trans.within_shuffled.CI.(multiple_trans.parameters.tpofinterest{i}) = [...
        prctile(multiple_trans.within_shuffled.(multiple_trans.parameters.tpofinterest{i}),2.5)... 
        prctile(multiple_trans.within_shuffled.(multiple_trans.parameters.tpofinterest{i}),97.5)];

end


%% Graph showing CI of history dependence and shuffled history dependence (within song)


for i = 1:length(multiple_trans.parameters.tpofinterest)
    if multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{i})(1) >...
            multiple_trans.within_shuffled.CI.(multiple_trans.parameters.tpofinterest{i})(2)
        multiple_trans.summary.is_bird_history_dependent.(multiple_trans.parameters.tpofinterest{i}) = 'yes';
        marker_color = [1 0 0];
    %the case if it is a linear song
    elseif multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{i})(1) == 1 && ...
            multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{i})(2) == 1
        multiple_trans.summary.is_bird_history_dependent.(multiple_trans.parameters.tpofinterest{i}) = 'yes';
        marker_color = [1 0 0];
    else
        multiple_trans.summary.is_bird_history_dependent.(multiple_trans.parameters.tpofinterest{i}) = 'no';
        marker_color = [0 0 0];
    end
    
    
    %figure for history dependence
    figure(), hold on
    
    %line for 95% CI of shuffled songs
    line([multiple_trans.within_shuffled.CI.(multiple_trans.parameters.tpofinterest{i})(1)...
        multiple_trans.within_shuffled.CI.(multiple_trans.parameters.tpofinterest{i})(2)],...
        [1 1],...
        'LineWidth', 2, 'Color', [.5 .5 .5])
    
    %line and marker for measured history dependence
    plot( multiple_trans.hist_dep.median.(multiple_trans.parameters.tpofinterest{i}),1.5,...
        'Marker', 'v', 'MarkerSize', 20, 'MarkerFaceColor', marker_color, 'MarkerEdgeColor', marker_color)
    line([multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{i})(1)...
        multiple_trans.hist_dep.CI.(multiple_trans.parameters.tpofinterest{i})(2)], [1.5 1.5],...
        'LineWidth', 2, 'Color', marker_color)
    set(gca,'YTick',[])
    ylim([0 2])
    title(['History dependence of ' multiple_trans.parameters.bird_name...
        ' for ' multiple_trans.parameters.tpofinterest{i} ' on ' multiple_trans.parameters.date_of_exp...
        '    (within song shuffle)'])
    xlabel(['History dependence (ex: | p(ab|ab) - p(ab|ac) - p(ab|ad) |)'])
    
    %saves figure
    saveas(figure(gcf), ['multi_sequence_' multiple_trans.parameters.date_of_analysis '_' multiple_trans.parameters.con_or_div '_'...
        multiple_trans.parameters.phrase ...
        '/Is_History_Dependent_(within_song_shuffle)_' multiple_trans.parameters.tpofinterest{i} '_' multiple_trans.parameters.date_of_exp], 'fig')
end

%% Display results and saves it all to a summary field

display(' ')
display(' ')
display('Transition Probability:  ')
display(multiple_trans.trans_prob.median)
display(' ')
display('Transition Probability Confidence Interval:  ')
display(multiple_trans.trans_prob.CI)

for i = 1:length(multiple_trans.parameters.tpofinterest)
    display(' ')
    display(['Transition probability for ' multiple_trans.parameters.tpofinterest{i} ' (songs kept intact):  '])
    display(multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{i}).median)
    display(['Standard deviation of transition probability for ' multiple_trans.parameters.tpofinterest{i} ' (songs kept intact):  '])
    display(multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{i}).sd)
%     display(['Coefficient of variation of transition probability for ' multiple_trans.parameters.tpofinterest{i} ' (songs kept intact):  '])
%     display(multiple_trans.trans_prob.rel_freq.(multiple_trans.parameters.tpofinterest{i}).cv)
end

display(' ')
display(' ')
display('Entropy: ')
display(multiple_trans.entropy.median)
display('Entropy CI: ')
display(multiple_trans.entropy.CI)

display(' ')
display(' ')
display('History dependence:  ')
display(multiple_trans.hist_dep.median)
display(' ')
display('History Dependence Confidence Interval:  ')
display(multiple_trans.hist_dep.CI)

display(' ')
display('Shuffled History Dependence Confidence Interval (within song shuffle): ')
display(multiple_trans.within_shuffled.CI)

display(' ')
display('Is history dependent?  ')
display(multiple_trans.summary.is_bird_history_dependent)

multiple_trans.summary.trans_prob.median = multiple_trans.trans_prob.median;
multiple_trans.summary.trans_prob.CI = multiple_trans.trans_prob.CI;

multiple_trans.summary.entropy.median = multiple_trans.entropy.median;
multiple_trans.summary.entropy.CI = multiple_trans.entropy.CI;

multiple_trans.summary.hist_dep.median = multiple_trans.hist_dep.median;
multiple_trans.summary.hist_dep.CI = multiple_trans.hist_dep.CI;

multiple_trans.summary.shuffled_hist_dep.CI = multiple_trans.within_shuffled.CI;

%% Making variable to look at multiple days when you are done for the day


if exist([multiple_trans.parameters.computer '/' multiple_trans.parameters.bird_name '/all_days_sequence_' multiple_trans.parameters.phrase], 'dir') ~= 7
    mkdir([multiple_trans.parameters.computer '/' multiple_trans.parameters.bird_name '/all_days_sequence_' multiple_trans.parameters.phrase])
else
end

save([multiple_trans.parameters.computer '/' multiple_trans.parameters.bird_name...
    '/all_days_sequence_' multiple_trans.parameters.phrase '/' multiple_trans.parameters.bird_name '_' ...
     multiple_trans.parameters.date_of_exp '_' multiple_trans.parameters.con_or_div], 'multiple_trans')



%% Save Data




% save(['multi_sequence_' multiple_trans.parameters.date_of_analysis '/multi_sequence_' multiple_trans.parameters.date_of_analysis])

% %move batch files to new folder, make sure the batch file name is not too
% %similar to the name of your bird
% if strcmpi(input(['Move batch files to multi_sequence_' multiple_trans.parameters.date_of_analysis ' folder? (y or n)  '], 's'),'y')
%     movefile([multiple_trans.parameters.batchfile(1:4) '*'],['multi_sequence_' multiple_trans.parameters.date_of_analysis])
% else
% end




    