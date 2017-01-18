%% LT 12/19/16 - ANOVA timecourse, shuffling transitions, etc.
function lt_neural_MultNeur_AnovaTcourse(NeuronDatabase)

% 1) multiple back (e.g. 2 back) -
% 2) positive control - ab vs. cd

window_fr = 0.02;
windshift_fr = 0.005;
min_trials_for_FR_and_ANOVA = 5; % for fr, actual data, and digrams shuffled.

wind_anova = 0.04; %
windshift_anova = 0.01;

max_gap_dur = 0.2; % for shuffled digrams and for actual data

%% ====


DATSTRUCT = struct;
NumNeurons = length(NeuronDatabase.neurons);


%% EXTRACT RAW DAT (syls and spikes)
for nn=1:NumNeurons
    disp(['neuron : ' num2str(nn)]);
    % 1) -- load dat
    cd(NeuronDatabase.global.basedir);
    
    % - find day folder
    dirdate=NeuronDatabase.neurons(nn).date;
    tmp=dir([dirdate '*']);
    assert(length(tmp)==1, 'PROBLEM - issue finding day folder');
    cd(tmp(1).name);
    
    % - load data for this neuron
    batchf=NeuronDatabase.neurons(nn).batchfile;
    channel_board=NeuronDatabase.neurons(nn).chan;
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);
    
    
    % ==== EXTRACT SPIKES, SYLS, FOR THIS NEURON
    DATSTRUCT.neuron(nn).song_labels = SongDat.AllLabels;
    DATSTRUCT.neuron(nn).song_onsets = SongDat.AllOnsets;
    DATSTRUCT.neuron(nn).song_offsets = SongDat.AllOffsets;
    
    DATSTRUCT.neuron(nn).spikes_cluster_class = NeurDat.spikes_cat.cluster_class;
    %     DATSTRUCT.neuron(nn).spikes_fs = NeurDat.metaDat(1).fs;
    
    DATSTRUCT.neuron(nn).params_batchf = Params.batchf;
    DATSTRUCT.neuron(nn).params_channel = Params.channel_board;
    
    DATSTRUCT.neuron(nn).metaDat = NeurDat.metaDat;
    DATSTRUCT.neuron(nn).metaDat = rmfield(DATSTRUCT.neuron(nn).metaDat, 'songDat');
end


%% CONVERGENT TRANSITION

for nn=1:NumNeurons
    
    allsyls = unique(DATSTRUCT.neuron(nn).song_labels);
    
    % --- format data for this neuron in manner useful for upcoming
    % analysesi.
    clear SongDat
    clear NeurDat
    clear Params
    SongDat.AllLabels = DATSTRUCT.neuron(nn).song_labels;
    SongDat.AllOffsets = DATSTRUCT.neuron(nn).song_offsets;
    SongDat.AllOnsets = DATSTRUCT.neuron(nn).song_onsets;
    NeurDat.spikes_cat.cluster_class = DATSTRUCT.neuron(nn).spikes_cluster_class;
    %         NeurDat.metaDat(1).fs = DATSTRUCT.neuron(nn).spikes_fs;
    NeurDat.metaDat = DATSTRUCT.neuron(nn).metaDat;
    
    % --- 1 figure for each neuron
    figcount=1;
    subplotrows=6;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    
    for i=1:length(allsyls)
        
        
        %% 1) ACTUAL CONVERGENT TRANS
        syl_master = allsyls(i);
        
        if strcmp(syl_master, '-')
            continue
        end
        
        % --- FIND ALL TRANSITIONS FOR THIS MASTER SYL
        inds = strfind(DATSTRUCT.neuron(nn).song_labels, syl_master);
        inds2 = inds(inds>1) - 1;
        
        syls_preceding = DATSTRUCT.neuron(nn).song_labels(inds2);
        syllist_slave = unique(syls_preceding);
        
        % ---------------------------------------------------------------------------------------------
        % === METHOD 1 used to be here - see below (1.1)
        % ==== METHOD 2
        % 1) get i) list of all transitions ii) sample size for each
        % (randomly drawn without replacement - i.e. if equals actual
        % sample size then will not be random)
        
        GroupLabels = {};
        GroupSampleSizes = []; % if use nan, then will take all. if use a number, then take random
        for ii=1:length(syllist_slave)
            syl_slave = syllist_slave(ii);
            
            if strcmp(syl_slave, '-')
                continue
            end
            
            regexpr_str = [syl_slave '(' syl_master ')'];
            
            GroupLabels = [GroupLabels regexpr_str];
            GroupSampleSizes = [GroupSampleSizes nan]; % to get all
        end
        
        % 2) EXTRACT ANOVA STRUCT
        % PREVIOUS version, before put into function, stored below (1.2)
        % added here: throws out transition if not enough trials
        predur = 0.15;
        postdur = 0.02;
        ANOVA_STRUCT = lt_neural_MultNeur_AnovaTcourse_sub1(GroupLabels, ...
            GroupSampleSizes, predur, postdur, SongDat, NeurDat, NeuronDatabase, nn, ...
            min_trials_for_FR_and_ANOVA);
        
        
        % ==== 1) PLOT, average fr time course for each category
        CategoriesList = unique({ANOVA_STRUCT.matchlabel});
        if length(CategoriesList)<2
            continue
        end
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['master syl ' syl_master]);
        fn_plot_fr(ANOVA_STRUCT, window_fr, windshift_fr, predur);
        
        
        % ==== 2) RUNNING ANOVA
        AllSpks_anova = {ANOVA_STRUCT.spk_Times};
        FactorLevels = {ANOVA_STRUCT.matchlabel};
        [Xall_ACTUAL, OmegaSquared_ACTUAL, ~, ~] = lt_neural_RunningAnova(AllSpks_anova, ...
            FactorLevels, wind_anova, windshift_anova);
        
        
        % ---- PLOT
        if (0) % skip becuase plotting again overlayed over shuffled ones.
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(['master syl ' syl_master '(ANOVA)']);
            
            plot(Xall_ACTUAL, OmegaSquared_ACTUAL, '-');
            axis tight
            line([predur predur], ylim, 'Color', 'k');
            ylim([-0.1 0.3]); lt_plot_zeroline
        end
        
        %% 2) FOR THIS SYL, GET POSITIVE CONTROL
        % i.e. replace each transition with a digram from the pool of all
        % digram types - match sample sizes. calculate anova. redo N times.
        
        % 1) === get list of all digrams present in his song.
        numsyls = length(SongDat.AllLabels);
        alldigrams = [SongDat.AllLabels(1:numsyls-1)', SongDat.AllLabels(2:numsyls)']; % get all digrams array
        alldigrams = mat2cell(alldigrams, ones(size(alldigrams,1), 1), 2); % one digram in each cell.
        
        % -- only keep instances in which gap in digram is < threshold dur
        digram_gaps = SongDat.AllOnsets(2:end)-SongDat.AllOffsets(1:end-1);
        assert(length(alldigrams) == length(digram_gaps), 'asdfasdf');
        inds = digram_gaps < max_gap_dur;
        alldigrams = alldigrams(inds); % done, kept those with short gaps;
        
        % -- get unique digrams
        alldigrams_unique = unique(alldigrams); % final, unique digrams
        
        Alldigrams_final = {}; % output
        Alldigrams_samplesize = []; % output
        for j=1:length(alldigrams_unique)
            digramcurr = alldigrams_unique{j};
            
            % -- remove all with dash
            if ~isempty(strfind(digramcurr, '-'))
                continue
            end
            
            % -- only keep those that occur at least N times
            if (0)
                % === VERSION 1, most accurate, but sample size is potentially larger than with
                % regexp, since the latter does not take overlapping (i.e.
                % bbbb, looking for bb, gives just 2 results). Either modify
                % code here to give regexp type output, or modify regexp to
                % work well for repeats. [this is only issue for repeats]
                % use lookaround? - e.g. 'b(?=b)' to get bb, but issue is can't
                % get it to work for getting tokens.
                samplesize = sum(strcmp(alldigrams, digramcurr));
            else
                % === VERSION 2 - same as regexp, so will not get larger sample
                % size than subsequent regexp code. use this until fix
                % regexp issue, described above.
                
                [startinds, ~, ~, ~]=regexp(SongDat.AllLabels, ...
                    digramcurr, 'start', 'end', 'match', 'tokenExtents');
                samplesize = length(startinds);
            end
            
            if  samplesize < min_trials_for_FR_and_ANOVA;
                continue
            end
            
            % -- if pass, then keep
            Alldigrams_final = [Alldigrams_final digramcurr];
            Alldigrams_samplesize = [Alldigrams_samplesize samplesize];
        end
        
        % 2 === for this convergent, choose random digram comparison
        % for actual data, distribtuion (across transitions) of sample
        % sizes
        UniqueLabels = unique({ANOVA_STRUCT.matchlabel});
        
        if (0)
            %         hfig = lt_figure; hold on;
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title('digrams chosen for shuffle');
            xlabel('shuff num'); ylabel('digram ind'); grid on;
            plot_digram_choices =1
        else
            plot_digram_choices = 0;
        end
        
        NCycles = 15; % cycles of getting surrogate data
        plotcols = lt_make_plot_colors(length(UniqueLabels), 0, 0);
        OmegaSquared_All= {};
        Xall_ALL={};
        for k=1:NCycles
            Digrams_Chosen = {};
            Digrams_Chosen_inds = [];
            Digrams_DesiredSampSize = [];
            for j=1:length(UniqueLabels)
                uniquelabel = UniqueLabels{j};
                samplesize = sum(strcmp({ANOVA_STRUCT.matchlabel}, uniquelabel));
                
                % --- for each real unique label, go through all the digrams and
                % pick out a random digram that has sample size at laest as large
                % as for this digram
                inds_potentialdigrams = find(Alldigrams_samplesize >= samplesize);
                
                tmp = 0;
                while tmp ==0
                    digram_chosen_ind =inds_potentialdigrams(randi(length(inds_potentialdigrams)));
                    % make sure this ind has not been chosen before
                    if any(Digrams_Chosen_inds == digram_chosen_ind)
                        tmp =0;
                        
                        if all(ismember(inds_potentialdigrams, Digrams_Chosen_inds))
                            % then potential digrams make a subset of the chosen digrams -
                            % i.e. no digram left to choose
                            
                            % OLD VERSION: - if there is no other unchosen digram with a
                            % large enough sample size, then downsample the
                            % actual data
                            %                             inds = find(strcmp({ANOVA_STRUCT.matchlabel}, uniquelabel));
                            %                             ind_to_remove = inds(randi(length(inds))); % pick a random ind to remove
                            %                             ANOVA_STRUCT(ind_to_remove) = [];
                            %
                            %                             samplesize = sum(strcmp({ANOVA_STRUCT.matchlabel}, uniquelabel)); % recalc sample size
                            %                             inds_potentialdigrams = find(Alldigrams_samplesize >= samplesize);
                            %                             disp(['samp size: '  num2str(samplesize)]);
                            
                            % NEW VERSION: downsampling is incorrect, since
                            % previous shuffles and actual data are using
                            % full sample. Simple skip this shuffle
                            % instance, and retry (note this down for
                            % diagnostics) - EXIT OUT OF THIS LOOP, AND
                            % RESTART FROM J=1
                            
                        end
                    else
                        tmp =1;
                    end
                end
                
                Digrams_Chosen = [Digrams_Chosen Alldigrams_final{digram_chosen_ind}];
                Digrams_Chosen_inds = [Digrams_Chosen_inds digram_chosen_ind];
                Digrams_DesiredSampSize = [Digrams_DesiredSampSize samplesize];
            end
            
            if plot_digram_choices == 1
                plot(k, Digrams_Chosen_inds, 'o');
            end
            
            % ==== EXTRACT ANOVA STRUCT
            GroupLabels = {};
            GroupSampleSizes = []; % if use nan, then will take all. if use a number, then take random
            for j=1:length(Digrams_Chosen)
                
                regexpr_str = [Digrams_Chosen{j}(1) '(' Digrams_Chosen{j}(2) ')'];
                
                GroupLabels = [GroupLabels regexpr_str];
                GroupSampleSizes = [GroupSampleSizes Digrams_DesiredSampSize(j)];
            end
            
            % 2) EXTRACT ANOVA STRUCT
            % PREVIOUS version, before put into function, stored below (1.2)
            % added here: throws out transition if not enough trials
            predur = 0.15;
            postdur = 0.02;
            ANOVA_STRUCT_SHUFF = lt_neural_MultNeur_AnovaTcourse_sub1(GroupLabels, ...
                GroupSampleSizes, predur, postdur, SongDat, NeurDat, NeuronDatabase, nn, ...
                min_trials_for_FR_and_ANOVA);
            
            % ===== ANOVA
            AllSpks_anova = {ANOVA_STRUCT_SHUFF.spk_Times};
            FactorLevels = {ANOVA_STRUCT_SHUFF.matchlabel};
            [Xall, OmegaSquared, ~, Pval] = lt_neural_RunningAnova(AllSpks_anova, ...
                FactorLevels, wind_anova, windshift_anova);
            
            OmegaSquared_All = [OmegaSquared_All OmegaSquared];
            Xall_ALL = [Xall_ALL Xall];
            
            
            %  ==== ACTIVATE TO PLOT FR AND ANOVA FOR EACH SHUFFLE REND
            if (0)
                % ===== FR DIFF
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                fn_plot_fr(ANOVA_STRUCT_SHUFF, window_fr, windshift_fr, predur)
                
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                plot(Xall, OmegaSquared, '-');
                axis tight
                line([predur predur], ylim, 'Color', 'k');
                ylim([-0.1 0.3]); lt_plot_zeroline
                pause
            end
            
            
        end
        
        %         pause;
        %         close(hfig);
        
        if ShuffleWentWell ==1
            % ================================== COMPARE ACTUAL OMEGA SQUARED
            % TIMECOURSE TO DISTRIBUTION OF SHUFFLED
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title('all shuffled anova');
            ylabel('omega^2'); xlabel('sec');
            % - 1) shuffled
            for j=1:length(OmegaSquared_All)
                plot(Xall_ALL{j}, OmegaSquared_All{j}, ':');
            end
            % - 2) actual
            plot(Xall_ACTUAL, OmegaSquared_ACTUAL, '-r');
            axis tight
            line([predur predur], ylim, 'Color', 'k');
            ylim([-0.1 0.5])
            lt_plot_zeroline;
            
            % - 3) p-val for actual
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title('prob of shuff > actual');
            
            xmax = min(cellfun(@length, Xall_ALL)); % convert cell to mat (only keeping inds that all datapoints share)
            tmp_Omega = nan(length(OmegaSquared_All), xmax);
            for j=1:length(OmegaSquared_All)
                tmp_Omega(j, :) = OmegaSquared_All{j}(1:xmax);
            end
            
            xmax = min([length(OmegaSquared_ACTUAL) xmax]); % actual vs. shuffled
            Pvals = [];
            for j=1:xmax
                tmp = sum(OmegaSquared_ACTUAL(j) < tmp_Omega(:, j))/size(tmp_Omega,1);
                Pvals = [Pvals tmp];
            end
            
            plot(Xall_ACTUAL(1:xmax), Pvals, 'g-o'); % pval
            axis tight; line([predur predur], ylim, 'Color', 'k');
            lt_plot_zeroline;
            ylim([0 1]);
        end
    end
    
    % ================= To do, collect all shuffled and actual anovas for
    % this neuron. Also do across all neurons.
    
end

end





function fn_plot_fr(ANOVA_STRUCT, window_fr, windshift_fr, predur)

CategoriesList= unique({ANOVA_STRUCT.matchlabel});

% -- calculate mean firing rates
%  first convert to cell of firing rates, do this once
%  for each category
plotcols = lt_make_plot_colors(length(CategoriesList), 0, 0);
for j=1:length(CategoriesList)
    currcat = CategoriesList{j};
    
    inds_trials = strcmp({ANOVA_STRUCT.matchlabel}, currcat);
    
    Yspks = {ANOVA_STRUCT(inds_trials).spk_Times};
    
    [xbin, ~, ~, ymean_hz, ysem_hz, ~, ystd_hz] = ...
        lt_neural_plotRastMean(Yspks, window_fr, windshift_fr, 0, '');
    
    shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{j}}, 1);
    %             plot(xbin, ymean_hz, '-')
    
    axis tight
    line([predur predur], ylim, 'Color', 'k');
end
end




%% SCRATCH

% === OLD METHOD 1 FOR CONVERGENT TRANSTIIONS (1.1)
%             % --- FIND ALL TRANSITIONS FOR THIS MASTER SYL
%             inds = strfind(DATSTRUCT.neuron(nn).song_labels, syl_master);
%             inds2 = inds(inds>1) - 1;
%
%             syls_preceding = DATSTRUCT.neuron(nn).song_labels(inds2);
%             syllist_slave = unique(syls_preceding);
%
%
%             clear ANOVA_STRUCT % holds, for this comparison, all data
%             % ----- Go thru all slaves and extract data for each one.
%             for ii=1:length(syllist_slave)
%                 syl_slave = syllist_slave(ii);
%
%                 if strcmp(syl_slave, '-')
%                     continue
%                 end
%
%                 regexpr_str = [syl_slave '(' syl_master ')'];
%                 predur = 0.15;
%                 postdur = 0.02;
%                 alignByOnset = 1;
%
%                 % - Extract data
%                 [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, '', ...
%                     regexpr_str, predur, postdur, alignByOnset, '');
%
%
%                 % ====== keep only spikes of a certain cluster number
%                 for iii = 1:length(SegmentsExtract)
%
%                     inds = SegmentsExtract(iii).spk_Clust == NeuronDatabase.neurons(nn).clustnum;
%                     SegmentsExtract(iii).spk_Clust = SegmentsExtract(iii).spk_Clust(inds);
%                     SegmentsExtract(iii).spk_Times = SegmentsExtract(iii).spk_Times(inds);
%                 end
%
%
%                 % ===== COLLECT ALL (for this master) INTO ONE STRUCTURE
%                 % important information for performing ANOVA - i) alignment
%                 % time, ii) category (for each trial)
%
%                 tmp_struct = SegmentsExtract; % only keep the essentials here
%                 tmp_struct = rmfield(tmp_struct, {'song_filename', 'song_datenum', ...
%                     'song_ind_in_batch', 'global_ontime_motifInclFlank', 'global_offtime_motifInclFlank', ...
%                     'fs', 'global_startind_motifonly', 'global_endind_motifonly', ...
%                     'global_tokenind_DatAlignedToOnsetOfThis', 'spk_Clust'});
%
%                 if ~exist('ANOVA_STRUCT', 'var')
%                     ANOVA_STRUCT = tmp_struct;
%                 else
%                     ANOVA_STRUCT = [ANOVA_STRUCT tmp_struct];
%                 end
%
%             end


% === (1.2)
%                 clear ANOVA_STRUCT % holds, for this comparison, all data
%                 for ii=1:length(GroupLabels)
%
%                     regexpr_str = GroupLabels{ii};
%                     predur = 0.15;
%                     postdur = 0.02;
%                     alignByOnset = 1;
%
%                     % - Extract data
%                     [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, '', ...
%                         regexpr_str, predur, postdur, alignByOnset, '');
%
%                     % ---- TAKE RANDOM SUBSET IF NEEDED
%                     if ~isnan(GroupSampleSizes(ii))
%                         % then take random subset equal to entry
%                         % IN PROGRESS
%                     end
%
%                     % ---- keep only spikes of a certain cluster number
%                     for iii = 1:length(SegmentsExtract)
%                         inds = SegmentsExtract(iii).spk_Clust == NeuronDatabase.neurons(nn).clustnum;
%                         SegmentsExtract(iii).spk_Clust = SegmentsExtract(iii).spk_Clust(inds);
%                         SegmentsExtract(iii).spk_Times = SegmentsExtract(iii).spk_Times(inds);
%                     end
%
%                     % ----- COLLECT ALL (for this master) INTO ONE STRUCTURE
%                     % important information for performing ANOVA - i) alignment
%                     % time, ii) category (for each trial)
%                     tmp_struct = SegmentsExtract; % only keep the essentials here
%                     tmp_struct = rmfield(tmp_struct, {'song_filename', 'song_datenum', ...
%                         'song_ind_in_batch', 'global_ontime_motifInclFlank', 'global_offtime_motifInclFlank', ...
%                         'fs', 'global_startind_motifonly', 'global_endind_motifonly', ...
%                         'global_tokenind_DatAlignedToOnsetOfThis', 'spk_Clust'});
%
%                     if ~exist('ANOVA_STRUCT', 'var')
%                         ANOVA_STRUCT = tmp_struct;
%                     else
%                         ANOVA_STRUCT = [ANOVA_STRUCT tmp_struct];
%                     end
%
%                 end
