function lt_neural_MultNeur_pseudopop(NeuronDatabase, TRANSMATRIX, branchtype, plotOnlySummary)


% minNeuronsPerBranch = 1; % only plot a branch point if N neurons have this labeled
minNperTran = 8; % min sample size
% motif_predur=0.13;
% motif_postdur=0.05;

spktimefield='spk_Times';
window = 0.02;
windshift = 0.005;
maxTransTime = 0.2; % in seconds, max time to allow for a transition - throws out data with longer trans.
%%
NumNeurons = length(TRANSMATRIX.neuron);


%% -- get list of all possible syls across all neurons
allsyls = [];
for i=1:NumNeurons
    allsyls = [allsyls TRANSMATRIX.neuron(i).sysInOrder];
end
allsyls = unique(allsyls);




%% ============================= COLLECT DAT FOR EACH BRANCH POINT.

DATSTRUCT = struct;

% for branchtype = {'conv', 'div'};
%     branchtype = branchtype{1};
for nn=1:NumNeurons
    disp(['neuron : ' num2str(nn)]);
    % 1) -- load dat
    try
        cd(NeuronDatabase.global.basedir);
    catch err
        cd(NeuronDatabase.neurons(nn).basedir);
    end
        
    dirdate=NeuronDatabase.neurons(nn).date;
    tmp=dir([dirdate '*']);
    if length(tmp)>1
        tmp = dir([dirdate '_' NeuronDatabase.neurons(nn).exptID '*']);
        assert(length(tmp) ==1,' daiosfhasiohfioawe');
    end
    cd(tmp(1).name);
    
%     % - find day folder
%     dirdate=NeuronDatabase.neurons(nn).date;
%     tmp=dir([dirdate '*']);
%     assert(length(tmp)==1, 'PROBLEM - issue finding day folder');
%     cd(tmp(1).name);
    
    % - load data for this neuron
    batchf=NeuronDatabase.neurons(nn).batchfile;
    channel_board=NeuronDatabase.neurons(nn).chan;
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);
    
    %     hsplots = cell(length(allsyls), 1);
    for i=1:length(allsyls)
        syl_master = allsyls(i); % this syl defines the brainch point - i..e if conv, then is syl2
        
        if strcmp(syl_master, '-')
            continue
        end
        
        Yspks_all={};
        SampSizes_all=[];
        % === for this neuron and branch point, get all transitions and
        % plot each one.
        for ii=1:length(allsyls)
            syl_slave = allsyls(ii);
            
            if strcmp(branchtype, 'conv')
                regexpr_str = [syl_slave '(' syl_master ')'];
                alignByOnset = 1;
                predur = 0.15;
                postdur = 0;
            elseif strcmp(branchtype, 'div')
                regexpr_str = ['(' syl_master ')' syl_slave];
                alignByOnset = 0;
                predur = 0.1;
                postdur = 0.1;
            end
            
            %             % ----- criteria to skip
            %             if strcmp(syl_slave, '-');
            %                 continue
            %             end
            % ------
            
            % ================ GET SMOOTHED FIRING RATE FOR THIS MOTIF
            % 2) -- extract dat
            disp(['extracting data for ' regexpr_str]);
            [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                regexpr_str, predur, postdur, alignByOnset, '');
            
            % -- throw out if not enough datapoints
            if isempty(fieldnames(SegmentsExtract))
                continue
            end
            
            % -- throw out if transition is too long
            tokeninds=[SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis];
            
            if strcmp(branchtype, 'conv');
                tmp = SongDat.AllOnsets(tokeninds) - SongDat.AllOffsets(tokeninds-1);
                inds_to_remove = tmp > maxTransTime;
            elseif strcmp(branchtype, 'div')
                tmp = SongDat.AllOnsets(tokeninds+1) - SongDat.AllOffsets(tokeninds);
                inds_to_remove = tmp > maxTransTime;
            end
            if sum(inds_to_remove) > 1
                disp(['Removed ' num2str(sum(inds_to_remove)) ' transtiions (too long)']);
            end
            SegmentsExtract(inds_to_remove) = [];
            
            
            
            % 3) ---------------------- get smoothed firing rate over all trials
            numtrials=length(SegmentsExtract);
            clustnum=NeuronDatabase.neurons(nn).clustnum;
            Yspks={};
            for k=1:numtrials
                inds=SegmentsExtract(k).spk_Clust==clustnum;
                spktimes=SegmentsExtract(k).(spktimefield)(inds);
                Yspks{k}=spktimes;
            end
            
            % -- convert to smoothed rate
            if all(cellfun(@isempty, Yspks))
                xbin = [];
                ymean_hz = [];
                ysem_hz = [];
                ystd_hz = [];
            else
                [xbin, ~, ~, ymean_hz, ysem_hz, ~, ystd_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
            end
            
            % ---- GET SYL time windows
            if strcmp(branchtype, 'conv')
                token_offtime = median(SongDat.AllOffsets([SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]) ...
                    - [SegmentsExtract.global_ontime_motifInclFlank]);
                %                 line([predur token_offtime], -rand(1)*10 - [5 5], 'Color', plotcol, 'LineWidth', 1.5)
                
            elseif strcmp(branchtype, 'div')
                token_ontime = median(SongDat.AllOnsets([SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]) ...
                    - [SegmentsExtract.global_ontime_motifInclFlank]);
                %                 line([token_ontime predur], -rand(1)*10 - [5 5], 'Color', plotcol, 'LineWidth', 1.5)
                
                
            end
            
            % ---- GET TIMING OF PRECEDING SYL
            globalind_presyl = [SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]-1;
            onsets_presyl = SongDat.AllOnsets(globalind_presyl);
            offsets_presyl = SongDat.AllOffsets(globalind_presyl);
            
            onsets_presyl = median(onsets_presyl - ...
                [SegmentsExtract.global_ontime_motifInclFlank]); % get relative to onset of extracted data
            offsets_presyl = median(offsets_presyl - ...
                [SegmentsExtract.global_ontime_motifInclFlank]);
            
            
            % ---- GET TIMING OF FOLLOWING SYL
            globalind_postsyl = [SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]+1;
            globalind_postsyl(globalind_postsyl>length(SongDat.AllOnsets))=[]; % remove things that don't have postsyl
            onsets_postsyl = SongDat.AllOnsets(globalind_postsyl);
            offsets_postsyl = SongDat.AllOffsets(globalind_postsyl);
            
            if isempty(onsets_postsyl)
                onsets_postsyl = nan;
                offsets_postsyl = nan;
            else
            onsets_postsyl = median(onsets_postsyl - ...
                [SegmentsExtract(1:length(onsets_postsyl)).global_ontime_motifInclFlank]);
            offsets_postsyl = median(offsets_postsyl - ...
                [SegmentsExtract(1:length(offsets_postsyl)).global_ontime_motifInclFlank]);
            end
            
            % =============== FOR THIS TRANSITION, SAVE MEAN FIRING CONTOUR
            DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).xbin = xbin; % i is master, ii is slave
            DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).ymean_hz = ymean_hz;
            DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).ysem_hz = ysem_hz;
            DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).ystd_hz = ystd_hz;
            DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N = length(SegmentsExtract);
            DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).Yspks_raw = Yspks;
            
            
            
            if strcmp(branchtype, 'conv')
                DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).token_onset_offset = [predur token_offtime];
            elseif strcmp(branchtype, 'div');
                DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).token_onset_offset = [token_ontime predur];
            end
            
            DATSTRUCT.neuron(nn).(branchtype).params.predur = predur;
            DATSTRUCT.neuron(nn).(branchtype).params.postdur = postdur;
            
            DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).presyl_onset_offset = ...
                [onsets_presyl offsets_presyl];
            DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).postsyl_onset_offset = ...
                [onsets_postsyl offsets_postsyl];
        end
        
        
        %         % ++++++++++++++++++++++++++ CONTROLS
        %         % ============ 1) for this branch, shuffle transitions - i.e.
        %         % assuming there is no difference in firing across transitions,
        %         % then what types of contours would we expect?
        %
        %         % Need:
        %         % raw spike train for each trial
        %         % sample sizes for each branch
        %         assert(sum(SampSizes_all) == length(Yspks_all));
        %         inds_shuff = randperm(length(Yspks_all));
        %
        %         DATSTRUCT.neuron(nn).(branchtype).transition_
        %
        %
        %
        %
        %
        %
        %         % ============ 2) take diff branches, how different is firing rate?
        %         % (i.e. max difference expected)
    end
end
% end

predur_global = predur;

%% ========== PLOT FIRIGN RATES AT ALL BRANCHES
%  NEED TO INCLUDE LINE FOR FLANKING SYLLABLES

if plotOnlySummary ==0
    
    for i=1:length(allsyls)
        syl_master = allsyls(i);
        if strcmp(syl_master, '-')
            continue
        end
        lt_figure; hold on;
        
        for nn=1:NumNeurons
            lt_subplot(8, 3, nn); hold on;
            title(['neuron ' num2str(nn)]);
            
            plotcols=lt_make_plot_colors(length(allsyls), 0, 0);
            for ii=1:length(allsyls)
                syl_slave = allsyls(ii);
                
                if strcmp(syl_slave, '-')
                    continue
                end
                
                %                 if size(DATSTRUCT.neuron(nn).(branchtype).transition, 1) >= i && ...
                %                         size(DATSTRUCT.neuron(nn).(branchtype).transition, 2) >= ii
                if DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N < minNperTran
                    continue
                end
                
                % ---- extract FR contour and plot
                xbin = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).xbin;
                ymean_hz = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).ymean_hz;
                ysem_hz = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).ysem_hz;
                token_onset_offset = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).token_onset_offset;
                presyl_onset_offset = ...
                    DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).presyl_onset_offset;
                postsyl_onset_offset = ...
                    DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).postsyl_onset_offset;
                
                if ~isempty(xbin)
                    shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{ii}}, 1);
                    lt_plot_text(xbin(end)+0.01, ymean_hz(end), syl_slave, plotcols{ii});
                    line(token_onset_offset, -rand(1)*10 - [5 5], 'Color', plotcols{ii}, 'LineWidth', 1.5)
                    %                     end
                    
                    % --- plot presyl and postsyl
                    if strcmp(branchtype, 'conv');
                        line(presyl_onset_offset, -rand(1)*10 - [5 5], 'Color', plotcols{ii}, 'LineWidth', 1.5)
                    elseif strcmp(branchtype, 'div');
                        line(presyl_onset_offset, -rand(1)*10 - [5 5], 'Color', plotcols{ii}, 'LineWidth', 1.5)
                        line(postsyl_onset_offset, -rand(1)*10 - [5 5], 'Color', plotcols{ii}, 'LineWidth', 1.5)
                    end
                    
                end
                
                axis tight
                
                % ---- Overlay vert line
                predur = DATSTRUCT.neuron(nn).(branchtype).params.predur;
                line([predur predur], ylim, 'Color', 'k');
            end
            
            
        end
            if strcmp(branchtype, 'conv')
                lt_subtitle(['CONV to ' syl_master]);
            elseif strcmp(branchtype, 'div');
                lt_subtitle(['DIV from ' syl_master]);
            end
    end
    if strcmp(input('close figs? (y or n)', 's'), 'y')
        close all;
    end
end

%% ======== FOR EACH BRANCH POINT, GET RUNNING AVERAGE OF FIRING RATE DEVIATIONS FOR ALL PAIRS OF TRANSITIONS
% =================== ACTUAL, GET ALL PAIRWISE D-PRIMES (BETWEEN ALL
% TRANSITIONS FOR A CERTAIN BRANCH POINT
if plotOnlySummary==0
    for i=1:length(allsyls)
        syl_master = allsyls(i);
        if strcmp(syl_master, '-')
            continue
        end
        lt_figure; hold on;
        
        FRdiffAcrossBranch={};
        XbinsAcrossBranch = {};
        for nn=1:NumNeurons
            lt_subplot(8, 3, nn); hold on;
            title(['neuron ' num2str(nn)]);
            ylabel('mean diff of FR');
            
            plotcols=lt_make_plot_colors(length(allsyls), 0, 0);
            AllSlaveContours = {};
            AllXbins= {};
            AllTokenOnOffsets = [];
            for ii=1:length(allsyls)
                syl_slave = allsyls(ii);
                if strcmp(syl_slave, '-')
                    continue
                end
                
                %             if size(DATSTRUCT.neuron(nn).(branchtype).transition, 1) >= i && ...
                %                     size(DATSTRUCT.neuron(nn).(branchtype).transition, 2) >= ii
                if DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N < minNperTran
                    continue
                end
                % ---- extract FR contour and plot
                xbin = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).xbin;
                ymean_hz = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).ymean_hz;
                ysem_hz = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).ysem_hz;
                token_onset_offset = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).token_onset_offset;
                
                AllXbins = [AllXbins; xbin];
                AllSlaveContours = [AllSlaveContours ymean_hz];
                AllTokenOnOffsets = [AllTokenOnOffsets; token_onset_offset];
                %             end
            end
            
            % ============= DO STUFF WITH THE SLAVE CONTOURS
            plot_on = 1; % plots for a given neuron, mean for this diff metric
            [FRdiffAcrossBranch, XbinsAcrossBranch] = fn_get_diff_from_slave_contour(...
                AllSlaveContours, plot_on, AllXbins, AllTokenOnOffsets, FRdiffAcrossBranch, XbinsAcrossBranch);
            
            % OLD VERSIONL:
            %         % ============= DO STUFF WITH THE SLAVE CONTOURS
            %         xmax = min(cellfun(@length, AllSlaveContours));
            %         % 1) get all paired differences of contour means
            %         numcontours = length(AllSlaveContours);
            %         if numcontours > 1
            %             tmp_alldiffs = [];
            %             for k=1:numcontours-1
            %                 for kk=k+1:numcontours
            %                     tmp_diff = abs(AllSlaveContours{k}(1:xmax) - AllSlaveContours{kk}(1:xmax));
            %                     tmp_alldiffs = [tmp_alldiffs; tmp_diff];
            %                 end
            %             end
            %
            %             % - PLOT
            %             plot(AllXbins{1}(1:xmax), mean(tmp_alldiffs, 1), '-o');
            %
            %             % ---- OVERLAY SYL ONSET OFFSET
            %             line(mean(AllTokenOnOffsets,1), -rand(1)*10 - [5 5], 'Color', 'k', 'LineWidth', 1.5)
            %             axis tight
            %
            %             % ---- Overlay vert line
            %             predur = DATSTRUCT.neuron(nn).(branchtype).params.predur;
            %             line([predur predur], ylim, 'Color', 'k');
            %
            %             lt_plot_zeroline;
            %
            %         end
        end
        
        if strcmp(branchtype, 'conv')
            lt_subtitle(['CONV to ' syl_master]);
        elseif strcmp(branchtype, 'div');
            lt_subtitle(['DIV from ' syl_master]);
        end
    end
end

%% ================= PLOT ALL DIFFS FOR EACH NEURON (I.E. EACH NEURON, ONE PLOT)
% COLLECT AND PLOT FOR EACH BRANCH POINT

if plotOnlySummary==0
    
    figcount=1;
    subplotrows=6;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    
    for nn=1:NumNeurons
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['neuron ' num2str(nn)]);
        ylabel('mean FR diff (hz)');
        xlabel('ms');
        
        FRdiffAcrossBranch = {};
        XbinsAcrossBranch = {};
        for i=1:length(allsyls)
            syl_master = allsyls(i);
            if strcmp(syl_master, '-')
                continue
            end
            
            % for this master (i.e. this branch) collect all contours, then
            % take difference between those contorus.
            AllSlaveContours = {};
            AllXbins= {};
            AllTokenOnOffsets = [];
            for ii=1:length(allsyls)
                syl_slave = allsyls(ii);
                if strcmp(syl_slave, '-')
                    continue
                end
                
                %             if size(DATSTRUCT.neuron(nn).(branchtype).transition, 1) >= i && ...
                %                     size(DATSTRUCT.neuron(nn).(branchtype).transition, 2) >= ii
                if DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N < minNperTran
                    continue
                end
                % ---- extract FR contour and plot
                xbin = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).xbin;
                ymean_hz = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).ymean_hz;
                %             plot(xbin, ymean_hz, '--');
                %                 ysem_hz = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).ysem_hz;
                token_onset_offset = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).token_onset_offset;
                
                AllXbins = [AllXbins xbin];
                AllSlaveContours = [AllSlaveContours ymean_hz];
                AllTokenOnOffsets = [AllTokenOnOffsets; token_onset_offset];
                %             end
            end
            
            % ============= DO STUFF WITH THE SLAVE CONTOURS
            plot_on = 1; % plots for a given neuron, mean for this diff metric
            [FRdiffAcrossBranch, XbinsAcrossBranch] = fn_get_diff_from_slave_contour(...
                AllSlaveContours, plot_on, AllXbins, AllTokenOnOffsets, FRdiffAcrossBranch, XbinsAcrossBranch);
            
            
            %         xmax = min(cellfun(@length, AllSlaveContours));
            %         % 1) get all paired differences of contour means
            %         numcontours = length(AllSlaveContours);
            %         if numcontours > 1
            %             tmp_alldiffs = [];
            %             for k=1:numcontours-1
            %                 for kk=k+1:numcontours
            %                     tmp_diff = abs(AllSlaveContours{k}(1:xmax) - AllSlaveContours{kk}(1:xmax));
            %                     tmp_alldiffs = [tmp_alldiffs; tmp_diff];
            %                 end
            %             end
            %
            %
            %             % ====== PLOT FOR THIS NEURON
            %             plot(AllXbins{1}(1:xmax), mean(tmp_alldiffs, 1), '-', 'Color', [0.2 0.7 0.7]); % this contains all diffs for this branch
            %
            %             % === PLOT LINE FOR SYL
            %             line(mean(AllTokenOnOffsets,1), -rand(1)*10 - [5 5], 'Color', [0.2 0.7 0.7], 'LineWidth', 1.5)
            %
            %             % ==== collect
            %             FRdiffAcrossBranch = [FRdiffAcrossBranch mean(tmp_alldiffs, 1)];
            %             XbinsAcrossBranch = [XbinsAcrossBranch AllXbins{1}];
            %         end
        end
        axis tight
        lt_plot_zeroline
        
        % ==== PLOT MEAN DIFF (across branches) for this neuron
        plotcol = [0.2 0.5 0.5];
        fn_plot_mean_diff_for_neuron(XbinsAcrossBranch, FRdiffAcrossBranch, plotcol)
        
        %     x2 = max(cellfun(@length, FRdiffAcrossBranch)); % not controlling for sample size
        %     tmptmp=nan(length(FRdiffAcrossBranch), x2);
        %     for lll = 1:length(FRdiffAcrossBranch)
        %         tmptmp(lll, 1:length(FRdiffAcrossBranch{lll})) = FRdiffAcrossBranch{lll};
        %     end
        %     xdiff = XbinsAcrossBranch{1}(2)-XbinsAcrossBranch{1}(1);
        %     xbins=XbinsAcrossBranch{1}(1):xdiff:(x2-1)*xdiff+XbinsAcrossBranch{1}(1);
        %     %     plot(xbins, nanmean(tmptmp,1), '-k');
        %     shadedErrorBar(xbins, nanmean(tmptmp, 1), lt_sem(tmptmp), {'-', 'Color', [0.2 0.5 0.5]},1 );
        
        % - plot, controlling for sample size (i.e. only overlapping regions)
        x2 = min(cellfun(@length, FRdiffAcrossBranch)); % controlling for sample size
        tmptmp=[];
        for lll = 1:length(FRdiffAcrossBranch)
            tmptmp = [tmptmp; FRdiffAcrossBranch{lll}(1:x2)];
        end
        shadedErrorBar(XbinsAcrossBranch{1}(1:x2), mean(tmptmp, 1), lt_sem(tmptmp), {'-k'}, 1);
        
        
        
        % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++===
        % ++++++++++++++++++++++++++++++++++++++++++++++ CONTROLS
        
        % ============ 1) for this branch, shuffle transitions - i.e.
        % assuming there is no difference in firing across transitions,
        % then what types of contours would we expect?
        FRdiffAcrossBranch = {};
        XbinsAcrossBranch = {};
        for i = 1:length(allsyls)
            disp(['master syl ' num2str(i)]);
            syl_master = allsyls(i);
            if strcmp(syl_master, '-')
                continue
            end
            
            % ------- COLLECT ALL RAW SPIKES (for transitions for this branch)
            Yspks_all = {};
            SampSizes_all = [];
            for ii=1:length(allsyls)
                syl_slave = allsyls(ii);
                if strcmp(syl_slave, '-')
                    continue
                end
                
                %             if size(DATSTRUCT.neuron(nn).(branchtype).transition, 1) >= i && ...
                %                     size(DATSTRUCT.neuron(nn).(branchtype).transition, 2) >= ii
                if DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N < minNperTran
                    continue
                end
                
                if isempty(DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).ymean_hz)
                    % use this as marker becuase these are the dartapoitns
                    % that went into actual data analysls.
                    continue
                end
                
                yspks_tmp = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).Yspks_raw;
                Yspks_all = [Yspks_all yspks_tmp];
                SampSizes_all = [SampSizes_all length(yspks_tmp)];
                %             end
            end
            assert(sum(SampSizes_all) == length(Yspks_all));
            
            % ----- FOR THIS BRANCH, SHUFFLE ALL THE SUBTRANSITIONS, AND GET
            % DIFFERENCE BETWEEN THEM
            inds_shuff = randperm(length(Yspks_all));
            numtrans = length(SampSizes_all);
            indlimits = [0 cumsum(SampSizes_all)];
            
            AllSlaveContours = {};
            AllXbins= {};
            for ii = 1:numtrans
                indstmp=inds_shuff(indlimits(ii)+1:indlimits(ii+1));
                Yspks_tmp = Yspks_all(indstmp);
                
                % -- calculate mean firing rates
                [xbin, ~, ~, ymean_hz, ysem_hz, ~, ystd_hz] = ...
                    lt_neural_plotRastMean(Yspks_tmp, window, windshift, 0, '');
                
                % -- save mean firing rate
                AllXbins = [AllXbins xbin];
                AllSlaveContours = [AllSlaveContours ymean_hz];
            end
            
            % ============= DO STUFF WITH THE SLAVE CONTOURS
            plot_on = 0; % plots for a given neuron, mean for this diff metric
            [FRdiffAcrossBranch, XbinsAcrossBranch] = fn_get_diff_from_slave_contour(...
                AllSlaveContours, plot_on, AllXbins, AllTokenOnOffsets, FRdiffAcrossBranch, XbinsAcrossBranch);
            
        end
        
        % ---- PLOT MEAN DIFF (across branches) for this neuron
        plotcol = [0.2 0.9 0.2];
        fn_plot_mean_diff_for_neuron(XbinsAcrossBranch, FRdiffAcrossBranch, plotcol)
        % ====================================================================
        
        
        % ============================================================================
        % ============= 2) take diff branches, how different is firing rate?
        % (i.e. max difference expected)
        % method: how diff are my data from case where all transitions to a branch
        % point are simply like different syllables? to test, for each branch
        % point, maybe has 5 transitions, each with sampel sizes [5 12 19 4 22].
        % Will get shuffled data by taking [5 12 19 4 22] from 5 different master
        % syllables (so get exact same sample size), each master collapsed across
        % all its slaves, randomly taken from the N syllables, with replacement.
        % Get difference of those contours..., then repeat for the next actual
        % branch.
        % MODIFICATION: instead of taking random syls, use the slave syl
        
        
        FRdiffAcrossBranch = {};
        XbinsAcrossBranch = {};
        for i = 1:length(allsyls)
            disp(['master syl ' num2str(i)]);
            syl_master = allsyls(i);
            if strcmp(syl_master, '-')
                continue
            end
            
            % slave data
            AllSlaveContours = {};
            AllXbins= {};
            %         AllTokenOnOffsets = [];
            for ii=1:length(allsyls)
                syl_slave = allsyls(ii);
                if strcmp(syl_slave, '-')
                    continue
                end
                
                if isempty(DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N)
                    continue
                end
                if DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N < minNperTran
                    continue
                end
                
                % ---- extract FR contour and plot
                sampsize = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N; % need to get this much of dat
                
                % ---- extract mean fr of slave (using slave as master)
                % - get all raw spike traces
                Yspks_all={};
                for iii=1:length(DATSTRUCT.neuron(nn).(branchtype).transition(ii,:))
                    Yspks_all = [Yspks_all ...
                        DATSTRUCT.neuron(nn).(branchtype).transition(ii,iii).Yspks_raw];
                end
                assert(length(Yspks_all) >= sampsize, 'not sure why, not enough sample of slave as new master');
                
                % - get random samplesize-sized inds
                inds_tmp = randperm(length(Yspks_all));
                inds_tmp = inds_tmp(1:sampsize);
                Yspks = Yspks_all(inds_tmp);
                
                % -- calculate mean firing rates
                [xbin, ~, ~, ymean_hz, ysem_hz, ~, ystd_hz] = ...
                    lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
                
                % -- output
                AllSlaveContours = [AllSlaveContours ymean_hz];
                AllXbins = [AllXbins xbin];
            end
            
            % ---- Do stuff with slave contours
            plot_on = 0;
            [FRdiffAcrossBranch, XbinsAcrossBranch] = fn_get_diff_from_slave_contour(...
                AllSlaveContours, plot_on, AllXbins, AllTokenOnOffsets, FRdiffAcrossBranch, XbinsAcrossBranch);
        end
        
        
        % ---- PLOT MEAN DIFF (across branches) for this neuron
        plotcol = 'r';
        fn_plot_mean_diff_for_neuron(XbinsAcrossBranch, FRdiffAcrossBranch, plotcol);
        % ====================================================================
        
    end
    
    lt_subtitle(['mean fr diff across branches (' branchtype ')']);
end

%% ========== TROUBLESHOOTING ABOVE - THIS IS LATEST - REPLACE ABOVE WITH THIS
% for a given neuron, and given master syl, plot 1) all actual transitions,
% 2) negative control - i.e. all shuffled, 3) positive control - i.e. all
% diff branches
close all;
NumShuffles = 1; % controls, to get distribution
deletefigs=0;

if plotOnlySummary == 0
if length(allsyls) * NumNeurons >40
    if strcmp(input('>40 figures, (syls * neurons), delete each batch of figs as they arise? (y or n)', 's'), 'y');
        deletefigs =1;
    end
end
end


syl_master_list = allsyls;
for sss = 1:length(syl_master_list);
    % syl_master_wanted = 'h';
    syl_master_wanted = syl_master_list(sss);
    
    if syl_master_wanted =='-';
        continue
    end
    
    
    for nn=1:NumNeurons
        lt_figure; hold on;
        hsplots = [];
        %     ylabel('mean FR diff (hz)');
        %     xlabel('ms');
        %
        %
        % ================================ 1) ALL ACTUAL TRANSITIONS
        FRdiffAcrossBranch = {};
        XbinsAcrossBranch = {};
        i= strfind(allsyls, syl_master_wanted);
        %     for i= [2 6];
        syl_master = allsyls(i);
        if strcmp(syl_master, '-')
            continue
        end
        hsplot = lt_subplot(3,3,1); hold on;
        title(['mean(SEM) fr (hz) - actual dat (neuron' ...
            num2str(nn) ' sylmaster' syl_master_wanted ')']);
        ylabel('actual transitions');
        xlabel('sec');
        line([predur_global predur_global], ylim, 'Color', 'k');

        hsplots = [hsplots hsplot];
        %         hsplotsy = [hsplotsy hsplot];
        % for this master (i.e. this branch) collect all contours, then
        % take difference between those contorus.
        AllSlaveContours = {};
        AllSlaveContours_std = {};
        AllXbins= {};
        AllTokenOnOffsets = [];
        AllPresylOnOff = [];
        AllPostsylOnOff = [];
        
        
        % for ANOVA analysis
        AllSpks_anova = {};
        AllSlaveLevels_anova = [];
        
        for ii=1:length(allsyls)
            syl_slave = allsyls(ii);
            if strcmp(syl_slave, '-')
                continue
            end
            
            %             if size(DATSTRUCT.neuron(nn).(branchtype).transition, 1) >= i && ...
            %                     size(DATSTRUCT.neuron(nn).(branchtype).transition, 2) >= ii
            if DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N < minNperTran
                continue
            end
            if isempty(DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N)
                continue
            end
            
            % ---- extract Spikes
            Yspks = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).Yspks_raw;
            % -- calculate mean firing rates
            [xbin, ~, ~, ymean_hz, ysem_hz, ~, ystd_hz] = ...
                lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
            
            % ---- extract FR contour and plot
            %             xbin = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).xbin;
            %             ymean_hz = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).ymean_hz;
            %             plot(xbin, ymean_hz, '--');
            %                 ysem_hz = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).ysem_hz;
            token_onset_offset = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).token_onset_offset;
            presyl_onset_offset = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).presyl_onset_offset;
            postsyl_onset_offset = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).postsyl_onset_offset;
            
            % --- PLOT
            plotcol = [rand rand rand];
            %             plot(xbin, ymean_hz, '-', 'Color', plotcol);
            shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcol}, 1);
            lt_plot_text(xbin(end), ymean_hz(end), allsyls(ii), plotcol);
            line(token_onset_offset,[-5+2*rand -5+2*rand], 'Color', plotcol)
                line([predur_global predur_global], ylim)
            
            try
                line(presyl_onset_offset, [-5+2*rand -5+2*rand], 'Color', plotcol);
            catch err
            end
            if strcmp(branchtype, 'div')
                line(postsyl_onset_offset, [-5+2*rand -5+2*rand], 'Color', plotcol);
            end
            
            axis tight
            lt_plot_zeroline;
            
            
            % --- SAVE
            AllXbins = [AllXbins xbin];
            AllSlaveContours = [AllSlaveContours ymean_hz];
            AllSlaveContours_std = [AllSlaveContours_std ystd_hz];
            AllTokenOnOffsets = [AllTokenOnOffsets; token_onset_offset];
            AllPresylOnOff = [AllPresylOnOff; presyl_onset_offset];
            AllPostsylOnOff = [AllPostsylOnOff; postsyl_onset_offset];
            
            
            % ----------------- 2) ANOVA EXTRACT
            AllSpks_anova = [AllSpks_anova Yspks];
            AllSlaveLevels_anova = [AllSlaveLevels_anova; ii*ones(length(Yspks),1)];
        end
        
        
        if length(AllSlaveContours)<2
            continue % go to next neruon
        end
        
        
        % ============= DO STUFF WITH THE SLAVE CONTOURS
        hsplot=lt_subplot(3,3,2); hold on;
        hsplots = [hsplots hsplot];
        diffmetric = 'dprime';
        plot_on = 1; % plots for a given neuron, mean for this diff metric
        %         diffmetric = 'frdiff';
        [FRdiffAcrossBranch, XbinsAcrossBranch, FRdiff_thisbranch, Xbin_forplot_thisbranch] = ...
            fn_get_diff_from_slave_contour(AllSlaveContours, plot_on, AllXbins, ...
            AllTokenOnOffsets, FRdiffAcrossBranch, XbinsAcrossBranch, diffmetric, AllSlaveContours_std);
        title(diffmetric);
        axis tight
        ylim([-1 2]);
        
        % ============= DO STUFF WITH RAW SPIKES
        % 1) ANOVA - each time bin, for a given master syl, get
        % variance explained by slave syl
        wind_anova = 0.04; %
        windshift_anova = 0.01;
        [Xall, OmegaSquared, EtaSquared, Pval] = lt_neural_RunningAnova(AllSpks_anova, ...
            AllSlaveLevels_anova, wind_anova, windshift_anova);
        % - plot
        hsplot = lt_subplot(3,3,3); hold on;
        title('anova');
        hsplots = [hsplots hsplot];
        plot(Xall, Pval, '-k');
        plot(Xall, OmegaSquared, '-r');
        lt_plot_annotation(1, 'bk: Pval; rd: Omega^2', 'k')
        axis tight
        lt_plot_zeroline
        ylim([-0.1 0.4]);
        line([predur_global predur_global], ylim, 'Color', 'k');

        % ================ SAVE CONTEXT METRICS
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).ActualDat.mastersyl = syl_master_wanted;
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).ActualDat.slave_syls = unique(AllSlaveLevels_anova);
        
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).ActualDat.(diffmetric).meanpairwise = FRdiff_thisbranch;
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).ActualDat.(diffmetric).xbins = Xbin_forplot_thisbranch;
        
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).ActualDat.anova.xbins = Xall;
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).ActualDat.anova.OmegaSquared = OmegaSquared;
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).ActualDat.anova.EtaSquared = EtaSquared;
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).ActualDat.anova.Pval = Pval;
        
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).ActualDat.AllTokenOnOffsets = AllTokenOnOffsets;
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).ActualDat.AllPresylOnOff = AllPresylOnOff;
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).ActualDat.AllPostsylOnOff = AllPostsylOnOff;
        
        
        
        
        % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++===
        % ++++++++++++++++++++++++++++++++++++++++++++++ CONTROLS
        
        % ============ 1) for this branch, shuffle transitions - i.e.
        % assuming there is no difference in firing across transitions,
        % then what types of contours would we expect?
        FRdiffAcrossBranch = {};
        XbinsAcrossBranch = {};
        i = strfind(allsyls, syl_master_wanted);
        disp(['master syl ' num2str(i)]);
        syl_master = allsyls(i);
        if strcmp(syl_master, '-')
            continue
        end
        
        %         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        %         title(['[shuff transitions (neg control): syl master: ' syl_master]);
        %         hsplots = [hsplots hsplot];
        %                 hsplotsy = [hsplotsy hsplot];
        
        % ------- COLLECT ALL RAW SPIKES (for transitions for this branch)
        Yspks_all = {};
        SampSizes_all = [];
        for ii=1:length(allsyls)
            syl_slave = allsyls(ii);
            if strcmp(syl_slave, '-')
                continue
            end
            
            %             if size(DATSTRUCT.neuron(nn).(branchtype).transition, 1) >= i && ...
            %                     size(DATSTRUCT.neuron(nn).(branchtype).transition, 2) >= ii
            if DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N < minNperTran
                continue
            end
            
            if isempty(DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N)
                continue
            end
            
            yspks_tmp = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).Yspks_raw;
            Yspks_all = [Yspks_all yspks_tmp];
            SampSizes_all = [SampSizes_all length(yspks_tmp)];
            %             end
        end
        assert(sum(SampSizes_all) == length(Yspks_all));
        
        
        % ----- FOR THIS BRANCH, SHUFFLE ALL THE SUBTRANSITIONS, AND GET
        % DIFFERENCE BETWEEN THEM
        tic
        FRdiff_shufflecollect = {};
        Pval_shufflecollect = {};
        OmegaSquared_shufflecollect = {};
        EtaSquared_shufflecollect = {};
        Xbin_frdiff= [];
        Xbin_anova = [];
        numtrans = length(SampSizes_all);
        indlimits = [0 cumsum(SampSizes_all)];
        
        for jjj = 1:NumShuffles
            inds_shuff = randperm(length(Yspks_all));
            
            AllSlaveContours = {};
            AllXbins= {};
            AllSlaveContours_std = {};
            % for ANOVA analysis
            AllSpks_anova = {};
            AllSlaveLevels_anova = [];
            
            for ii = 1:numtrans
                indstmp=inds_shuff(indlimits(ii)+1:indlimits(ii+1));
                Yspks_tmp = Yspks_all(indstmp);
                
                % -- calculate mean firing rates
                [xbin, ~, ~, ymean_hz, ysem_hz, ~, ystd_hz] = ...
                    lt_neural_plotRastMean(Yspks_tmp, window, windshift, 0, '');
                
                % == plot
                if jjj==1 % only plot one example
                    hsplot = lt_subplot(3,3,4); hold on;
                    hsplots = [hsplots hsplot];
                    ylabel('context-shuffled');
                    plotcol = [rand rand rand];
                    %             plot(xbin, ymean_hz, '-', 'Color', plotcol);
                    shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcol}, 1);
                    axis tight
                    lt_plot_zeroline
                end
                
                % -- save mean firing rate
                AllXbins = [AllXbins xbin];
                AllSlaveContours = [AllSlaveContours ymean_hz];
                AllSlaveContours_std = [AllSlaveContours_std ystd_hz];
                % ----------------- 2) ANOVA EXTRACT
                AllSpks_anova = [AllSpks_anova Yspks_tmp];
                AllSlaveLevels_anova = [AllSlaveLevels_anova; ii*ones(length(Yspks_tmp),1)];
            end
            
            
            % ============= DO STUFF WITH THE SLAVE CONTOURS
            plot_on = 0; % plots for a given neuron, mean for this diff metric
            diffmetric = 'dprime';
            %         diffmetric = 'frdiff';
            [~, ~, FRdiff_thisbranch, Xbin_forplot_thisbranch] ...
                = fn_get_diff_from_slave_contour(AllSlaveContours, plot_on, AllXbins, ...
                AllTokenOnOffsets, FRdiffAcrossBranch, XbinsAcrossBranch, diffmetric, AllSlaveContours_std);
            
            FRdiff_shufflecollect = [FRdiff_shufflecollect FRdiff_thisbranch];
            if length(Xbin_forplot_thisbranch)>length(Xbin_frdiff)
                Xbin_frdiff = Xbin_forplot_thisbranch;
            end
            
            % ============= DO STUFF WITH RAW SPIKES
            % 1) ANOVA - each time bin, for a given master syl, get
            % variance explained by slave syl
            [Xall, OmegaSquared, EtaSquared, Pval] = lt_neural_RunningAnova(AllSpks_anova, ...
                AllSlaveLevels_anova, wind_anova, windshift_anova);
            
            Pval_shufflecollect = [Pval_shufflecollect Pval];
            OmegaSquared_shufflecollect = [OmegaSquared_shufflecollect OmegaSquared];
            EtaSquared_shufflecollect = [EtaSquared_shufflecollect EtaSquared];
            if length(Xall) > length(Xbin_anova)
                Xbin_anova = Xall;
            end
            
        end
        toc
        line([predur_global predur_global], ylim, 'Color', 'k');

        
        % ======================= PLOT:
        hsplot =  lt_subplot(3,3,5); hold on;
        hsplots = [hsplots hsplot];
        yall=[];
        xmax = min(cellfun(@length, FRdiff_shufflecollect));
        for jjj=1:length(FRdiff_shufflecollect)
            xx = length(FRdiff_shufflecollect{jjj});
            plot(Xbin_frdiff(1:xx)', FRdiff_shufflecollect{jjj}, '-', 'Color', [0.7 0.7 0.7]);
            yall = [yall; FRdiff_shufflecollect{jjj}(1:xmax)];
        end
        %          errbar = prctile(yall(:,1:xmax), [97.5 2.5]);
%         shadedErrorBar(Xbin_frdiff(1:xmax), mean(yall(:,1:xmax),1), ...
%             2*std(yall, 0, 1), {'Color', 'b'},1);
        axis tight
        ylim([-1 2]);
        lt_plot_zeroline;
        line([predur_global predur_global], ylim, 'Color', 'k');

        % anova
        hsplot = lt_subplot(3,3,6); hold on;
        hsplots = [hsplots hsplot];
        xmax = min(cellfun(@length, Pval_shufflecollect));
        yall_p=[];
        yall_omega=[];
        for jjj=1:length(Pval_shufflecollect)
            xx = length(Pval_shufflecollect{jjj});
            plot(Xbin_anova(1:xx), Pval_shufflecollect{jjj}, '-', 'Color', [0.7 0.7 0.7]);
            plot(Xbin_anova(1:xx), OmegaSquared_shufflecollect{jjj}, '-', 'Color', [1 0.7 0.7]);
            
            yall_p = [yall_p; Pval_shufflecollect{jjj}(1:xmax)];
            yall_omega = [yall_omega; OmegaSquared_shufflecollect{jjj}(1:xmax)];
            
        end
%         shadedErrorBar(Xbin_anova(1:xmax), mean(yall_p, 1), 2*std(yall_p, 0,1), ...
%             {'Color', 'k'}, 1);
%         shadedErrorBar(Xbin_anova(1:xmax), mean(yall_omega, 1), 2*std(yall_omega, 0,1), ...
%             {'Color', 'r'}, 1);
        
        lt_plot_annotation(1, 'bk: Pval; rd: Omega^2', 'k')
        ylim([0 1]);
        axis tight
        ylim([-0.1 0.4]);
        line([predur_global predur_global], ylim, 'Color', 'k');

        % ================ SAVE CONTEXT METRICS
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).Control_ContextShuff.(diffmetric).meanpairwise_MultShuff = FRdiff_shufflecollect;
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).Control_ContextShuff.(diffmetric).xbins_MultShuff = Xbin_frdiff;
        
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).Control_ContextShuff.anova.xbins_MultShuff = Xbin_anova;
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).Control_ContextShuff.anova.OmegaSquared_MultShuff = OmegaSquared_shufflecollect;
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).Control_ContextShuff.anova.EtaSquared_shufflecollect = EtaSquared_shufflecollect;
        
        %                 DATSTRUCT.neuron(nn).(branchtype).branch_master(i).Control_ContextShuff.anova.EtaSquared = EtaSquared;
        DATSTRUCT.neuron(nn).(branchtype).branch_master(i).Control_ContextShuff.anova.Pval_MultShuff = Pval_shufflecollect;
        
        
        
        % ============================================================================
        % ============= 2) take diff branches, how different is firing rate?
        % (i.e. max difference expected)
        % method: how diff are my data from case where all transitions to a branch
        % point are simply like different syllables? to test, for each branch
        % point, maybe has 5 transitions, each with sampel sizes [5 12 19 4 22].
        % Will get shuffled data by taking [5 12 19 4 22] from 5 different master
        % syllables (so get exact same sample size), each master collapsed across
        % all its slaves, randomly taken from the N syllables, with replacement.
        % Get difference of those contours..., then repeat for the next actual
        % branch.
        % MODIFICATION (done): instead of taking random syls, use the slave syl
        
        FRdiffAcrossBranch = {};
        XbinsAcrossBranch = {};
        i = strfind(allsyls, syl_master_wanted);
        
        syl_master = allsyls(i);
        if strcmp(syl_master, '-')
            continue
        end
        
        %         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        %         title('diff syls control');
        %         hsplots = [hsplots hsplot];
        %                 hsplotsy = [hsplotsy hsplot];
        %
        %
        
        % ===========================
        % GET RANDOM SUBSET OF DATA MULTIPEL TIMES
        FRdiff_shufflecollect = {};
        Pval_shufflecollect = {};
        OmegaSquared_shufflecollect = {};
        EtaSquared_shufflecollect = {};
        Xbin_frdiff= [];
        Xbin_anova = [];
        
        for jjj=1:NumShuffles
            
            % slave data
            AllSlaveContours = {};
            AllXbins= {};
            AllSlaveContours_std = {};
            
            % for ANOVA analysis
            AllSpks_anova = {};
            AllSlaveLevels_anova = [];
            
            
            for ii=1:length(allsyls)
                syl_slave = allsyls(ii);
                if strcmp(syl_slave, '-')
                    continue
                end
                
                if isempty(DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N)
                    continue
                end
                if DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N < minNperTran
                    continue
                end
                
                % ---- extract FR contour and plot
                sampsize = DATSTRUCT.neuron(nn).(branchtype).transition(i, ii).N; % need to get this much of dat
                
                % ---- extract mean fr of slave (using slave as master)
                % - get all raw spike traces
                Yspks_all={};
                for iii=1:length(DATSTRUCT.neuron(nn).(branchtype).transition(ii,:))
                    
                    
                    Yspks_all = [Yspks_all ...
                        DATSTRUCT.neuron(nn).(branchtype).transition(ii,iii).Yspks_raw];
                end
                assert(length(Yspks_all) >= sampsize, 'not sure why, not enough sample of slave as new master');
                
                % - get random samplesize-sized inds
                inds_tmp = randperm(length(Yspks_all));
                inds_tmp = inds_tmp(1:sampsize);
                Yspks = Yspks_all(inds_tmp);
                
                % -- calculate mean firing rates
                [xbin, ~, ~, ymean_hz, ysem_hz, ~, ystd_hz] = ...
                    lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
                
                % === plot
                if jjj==1 % only plot one example
                    hsplot = lt_subplot(3,3, 7); hold on;
                    ylabel('slave syls (pos control)');
                    hsplots = [hsplots hsplot];
                    plotcol = [rand rand rand];
                    %             plot(xbin, ymean_hz, '-', 'Color', plotcol);
                    shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcol}, 1);
                    lt_plot_text(xbin(end), ymean_hz(end), allsyls(ii), plotcol);
                    axis tight
                    lt_plot_zeroline
                end
                
                % -- output
                AllSlaveContours = [AllSlaveContours ymean_hz];
                AllXbins = [AllXbins xbin];
                AllSlaveContours_std = [AllSlaveContours_std ystd_hz];
                
                
                % ----------------- 2) ANOVA EXTRACT
                AllSpks_anova = [AllSpks_anova Yspks];
                AllSlaveLevels_anova = [AllSlaveLevels_anova; ii*ones(length(Yspks),1)];
                
            end
            
            % ---- Do stuff with slave contours
            plot_on = 0;
            diffmetric = 'dprime';
            %         diffmetric = 'frdiff';
            [~, ~, FRdiff_thisbranch, Xbin_forplot_thisbranch] = fn_get_diff_from_slave_contour(...
                AllSlaveContours, plot_on, AllXbins, AllTokenOnOffsets, FRdiffAcrossBranch, ...
                XbinsAcrossBranch, diffmetric, AllSlaveContours_std);
            
            FRdiff_shufflecollect = [FRdiff_shufflecollect FRdiff_thisbranch];
            if length(Xbin_forplot_thisbranch)>length(Xbin_frdiff)
                Xbin_frdiff = Xbin_forplot_thisbranch;
            end
            
            
            % ============= DO STUFF WITH RAW SPIKES
            % 1) ANOVA - each time bin, for a given master syl, get
            % variance explained by slave syl
            [Xall, OmegaSquared, EtaSquared, Pval] = lt_neural_RunningAnova(AllSpks_anova, ...
                AllSlaveLevels_anova, wind_anova, windshift_anova);
            
            Pval_shufflecollect = [Pval_shufflecollect Pval];
            OmegaSquared_shufflecollect = [OmegaSquared_shufflecollect OmegaSquared];
            EtaSquared_shufflecollect = [EtaSquared_shufflecollect EtaSquared];
            if length(Xall) > length(Xbin_anova)
                Xbin_anova = Xall;
            end
            
        end
        
        % =========================
        % - plot - fr diff (e.g. dprime)
        hsplot = lt_subplot(3,3,8); hold on;
        hsplots = [hsplots hsplot];
        yall=[];
        xmax = min(cellfun(@length, FRdiff_shufflecollect));
        for jjj=1:length(FRdiff_shufflecollect)
            xx = length(FRdiff_shufflecollect{jjj});
            plot(Xbin_frdiff(1:xx), FRdiff_shufflecollect{jjj}, '-', 'Color', [0.7 0.7 0.7]);
            yall = [yall; FRdiff_shufflecollect{jjj}(1:xmax)];
        end
        %          errbar = prctile(yall(:,1:xmax), [97.5 2.5]);
%         shadedErrorBar(Xbin_frdiff(1:xmax), mean(yall(:,1:xmax),1), ...
%             2*std(yall, 0, 1), {'Color', 'b'},1);
        axis tight
        ylim([-1 2]);
        lt_plot_zeroline;
        line([predur_global predur_global], ylim, 'Color', 'k');

        % - anova
        hsplot = lt_subplot(3,3,9); hold on;
        hsplots = [hsplots hsplot];
        title('anova');
        xmax = min(cellfun(@length, Pval_shufflecollect));
        yall_p=[];
        yall_omega=[];
        for jjj=1:length(Pval_shufflecollect)
            xx = length(Pval_shufflecollect{jjj});
            plot(Xbin_anova(1:xx), Pval_shufflecollect{jjj}, '-', 'Color', [0.7 0.7 0.7]);
            plot(Xbin_anova(1:xx), OmegaSquared_shufflecollect{jjj}, '-', 'Color', [1 0.7 0.7]);
            
            yall_p = [yall_p; Pval_shufflecollect{jjj}(1:xmax)];
            yall_omega = [yall_omega; OmegaSquared_shufflecollect{jjj}(1:xmax)];
            
        end
%         shadedErrorBar(Xbin_anova(1:xmax), mean(yall_p, 1), 2*std(yall_p, 0,1), ...
%             {'Color', 'k'}, 1);
%         shadedErrorBar(Xbin_anova(1:xmax), mean(yall_omega, 1), 2*std(yall_omega, 0,1), ...
%             {'Color', 'r'}, 1);
        
        lt_plot_annotation(1, 'bk: Pval; rd: Omega^2', 'k')
        ylim([0 1]);
        axis tight
        ylim([-0.1 0.4]);
        line([predur_global predur_global], ylim, 'Color', 'k');

        
        % ================ SAVE CONTEXT METRICS
        DATSTRUCT.neuron(nn).(branchtype).branch_master...
            (i).ControlPos_SlaveSyls.(diffmetric).meanpairwise_MultShuff = FRdiff_shufflecollect;
        DATSTRUCT.neuron(nn).(branchtype).branch_master...
            (i).ControlPos_SlaveSyls.(diffmetric).xbins_MultShuff = Xbin_frdiff;
        
        DATSTRUCT.neuron(nn).(branchtype).branch_master...
            (i).ControlPos_SlaveSyls.anova.xbins_MultShuff = Xbin_anova;
        DATSTRUCT.neuron(nn).(branchtype).branch_master...
            (i).ControlPos_SlaveSyls.anova.OmegaSquared_MultShuff = OmegaSquared_shufflecollect;
        DATSTRUCT.neuron(nn).(branchtype).branch_master...
            (i).ControlPos_SlaveSyls.anova.EtaSquared_shufflecollect = EtaSquared_shufflecollect;
        %                 DATSTRUCT.neuron(nn).(branchtype).branch_master(i).Control_ContextShuff.anova.EtaSquared = EtaSquared;
        DATSTRUCT.neuron(nn).(branchtype).branch_master...
            (i).ControlPos_SlaveSyls.anova.Pval_MultShuff = Pval_shufflecollect;
        
        
        
        
        
        linkaxes(hsplots, 'x');
        % linkaxes(hsplotsy, 'y');
        lt_subtitle(['neuron ' num2str(nn) '; syl master: ' syl_master_wanted]);
        
    end
    
    if plotOnlySummary ==1 % then close automatically
        close all;
    elseif deletefigs==1
        disp('PRESS ANYTHING TO DELETE FIG AND CONTINUE');
        pause
        close all
    end
end



%% ====== JUST ANOVA R2 [EACH SUBPLOT, ONE NEURON]
% TO DO: get empirical p-val [compared to shuffle distribution] 
% at each time point and plot that
alpha = 0.05; % for anova p-val

figcount=1;
subplotrows=3;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];

plotActualDat = 1;

for nn=1:NumNeurons
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['neuron ' num2str(nn)]);
    ylabel('eta^2');
    
    numBranches = length(DATSTRUCT.neuron(nn).(branchtype).branch_master);
    plotcols = lt_make_plot_colors(numBranches, 0, 0);
    
    % --- get median presyl and postsyl timing (for this neuron)
    Presyl_onoff = [];
    for bbb = 1:length(DATSTRUCT.neuron(nn).(branchtype).branch_master);
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(bbb).ActualDat)
            continue
        end
        Presyl_onoff = [Presyl_onoff; ...
            DATSTRUCT.neuron(nn).(branchtype).branch_master(bbb).ActualDat.AllPresylOnOff];
    end
    
    if plotActualDat==1
    % -- 1) ACTUAL DATA
    plotcol = 'k';
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.xbins;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.EtaSquared;
        plot(xbins, yval, '-', 'Color', plotcol, 'LineWidth', 2);
    end
    
    else
        % -- 3) POS CONTROL
    plotcol = 'b';
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.EtaSquared_shufflecollect;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, yval, '-', 'Color', plotcol, 'LineWidth', 2);
    end
    end
    
    
    % -- 3) NEG CONTROL
    plotcol = 'r';
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.EtaSquared_shufflecollect;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, yval, '-', 'Color', plotcol, 'LineWidth', 1);
    end
%     axis tight
%     ylim([-0.1 0.4])
%     lt_plot_zeroline
%     line([predur_global predur_global], ylim);
%     
    
    
    axis tight
    ylim([-0.1 0.4])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
    xmean = mean(Presyl_onoff, 1);
    line(xmean, [-0.05 -0.05], 'Color', [0.3 0.3 0.8], 'LineWidth', 2);
    
    
end

%% ====== JUST ANOVA R2 [EACH SUBPLOT, ONE BRANCH]
% TO DO: get empirical p-val [compared to shuffle distribution] 
% at each time point and plot that
alpha = 0.05; % for anova p-val

figcount=1;
subplotrows=3;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];

plotActualDat = 1; % NOTE: CHANGE THIS TO CHOOSE WHETHER PLOT ACTUAL DAT OR POS CONTROL!!!!!

numBranches = length(allsyls);
for j=1:numBranches
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['branch ' allsyls(j)]);
    ylabel('eta^2');
    
    % --- get median presyl and postsyl timing (for this neuron)
    %     Presyl_onoff = [];
    %     for bbb = 1:length(DATSTRUCT.neuron(1).(branchtype).branch_master);
    %         if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(bbb).ActualDat)
    %             continue
    %         end
    %         Presyl_onoff = [Presyl_onoff; ...
    %             DATSTRUCT.neuron(nn).(branchtype).branch_master(bbb).ActualDat.AllPresylOnOff];
    %     end
    
    if plotActualDat==1
        % -- 1) ACTUAL DATA
        plotcol = 'k';
        for nn=1:NumNeurons
            if length(DATSTRUCT.neuron(nn).(branchtype).branch_master)<j
                continue
            end
            if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
                continue
            end
            
            xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.xbins;
            yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.EtaSquared;
            plot(xbins, yval, '-', 'Color', plotcol, 'LineWidth', 2);
        end
        
    else
        % -- 3) POS CONTROL
        plotcol = 'b';
        for nn=1:NumNeurons
            if length(DATSTRUCT.neuron(nn).(branchtype).branch_master)<j
                continue
            end
            if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
                continue
            end
            
            xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.xbins_MultShuff;
            yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.EtaSquared_shufflecollect;
            xmax = min(cellfun(@length, yval_cell));
            yval = [];
            for jj=1:length(yval_cell)
                yval = [yval; yval_cell{jj}(1:xmax)];
            end
            xbins = xbins(1:xmax);
            yval = mean(yval, 1);
            plot(xbins, yval, '-', 'Color', plotcol, 'LineWidth', 2);
        end
    end
    
    
    % -- 3) NEG CONTROL
    plotcol = 'r';
    for nn=1:NumNeurons
        if length(DATSTRUCT.neuron(nn).(branchtype).branch_master)<j
            continue
        end
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.EtaSquared_shufflecollect;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, yval, '-', 'Color', plotcol, 'LineWidth', 1);
    end
    %     axis tight
    %     ylim([-0.1 0.4])
    %     lt_plot_zeroline
    %     line([predur_global predur_global], ylim);
    %
    
    
    axis tight
    ylim([-0.1 0.4])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
    xmean = mean(Presyl_onoff, 1);
    line(xmean, [-0.05 -0.05], 'Color', [0.3 0.3 0.8], 'LineWidth', 2);
end

%% ====== COMBINE ALL THOSE CONTEXT METRICS ACROSS ALL BRANCHES FOR A GIVEN NEURON [FIG = NEURON]

alpha = 0.05; % for anova p-val

for nn=1:NumNeurons
    figcount=1;
    subplotrows=3;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    
    numBranches = length(DATSTRUCT.neuron(nn).(branchtype).branch_master);
    
    plotcols = lt_make_plot_colors(numBranches, 0, 0);
    
    % --- get median presyl and postsyl timing (for this neuron)
    Presyl_onoff = [];
    for bbb = 1:length(DATSTRUCT.neuron(nn).(branchtype).branch_master);
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(bbb).ActualDat)
            continue
        end         
        Presyl_onoff = [Presyl_onoff; ...
            DATSTRUCT.neuron(nn).(branchtype).branch_master(bbb).ActualDat.AllPresylOnOff];
    end
    
    
    
    % === ACTUAL DATA
    % -- 1) d-prime
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['dprime, NEURON ' num2str(nn)]);
    ylabel('ACTUAL dat');
    xlabel('sec');
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.dprime.xbins;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.dprime.meanpairwise;
        plot(xbins, yval, '-', 'Color', plotcols{j});
    end
    axis tight
    ylim([-1 2])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
    % -- line for presyl
    xstd = std(Presyl_onoff, 0, 1);
    xmean = mean(Presyl_onoff, 1);
    line(xmean, [-0.1 -0.1], 'Color', [0.3 0.3 0.8], 'LineWidth', 2);
%     line([xmean(1)-xstd(1) xmean(1)+xstd(1)], [-0.1 -0.1], 'LineWidth', 2);
%     plot(Presyl_onoff(:,1), -0.5, 'o');

        
    % -- 2) anova pval
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('log10(p-val)');
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.xbins;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.Pval;
        plot(xbins, log10(yval), '-', 'Color', plotcols{j});
    end
    axis tight
    ylim([-5 0]);
    line([predur_global predur_global], ylim);
    line(xlim, log10([alpha alpha]), 'Color', 'r', 'LineStyle', '--');
        xmean = mean(Presyl_onoff, 1);
    line(xmean, [-4 -4], 'Color', [0.3 0.3 0.8], 'LineWidth', 2);

    % proportion branches with signficant p-val
    
    % -- 3) anova - r2
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('eta^2');
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.xbins;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.EtaSquared;
        plot(xbins, yval, '-', 'Color', plotcols{j});
    end
    axis tight
    ylim([-0.1 0.4])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
        xmean = mean(Presyl_onoff, 1);
    line(xmean, [-0.05 -0.05], 'Color', [0.3 0.3 0.8], 'LineWidth', 2);

    
    % ============= NEG CONTROL (CONTEXT SHUFFLED)
    % -- 4) d-prime
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('mean x-shuff');
    ylabel('NEG CTRL (CTXT SHUFFLE)');
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.dprime.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.dprime.meanpairwise_MultShuff;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, yval, '-', 'Color', plotcols{j});
    end
    axis tight
    ylim([-1 2])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
    
    
    % -- 2) anova pval
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('mean x-shuff');
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.Pval_MultShuff;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, log10(yval), '-', 'Color', plotcols{j});
    end
    axis tight
    ylim([-5 0]);
    line([predur_global predur_global], ylim);
    line(xlim, log10([alpha alpha]), 'Color', 'r', 'LineStyle', '--');
    
    % proportion branches with signficant p-val
    
    % -- 3) anova - r2
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.EtaSquared_shufflecollect;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, yval, '-', 'Color', plotcols{j});
    end
    axis tight
    ylim([-0.1 0.4])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
    
    
    % ================= POS CONTROL (DIFF SYLS)
    % -- 4) d-prime
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('mean x-shuff');
    ylabel('POST CTRL (SLAVE SYls)');
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.dprime.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.dprime.meanpairwise_MultShuff;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, yval, '-', 'Color', plotcols{j});
    end
    axis tight
    ylim([-1 2])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
    
    
    % -- 2) anova pval
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('mean x-shuff');
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.Pval_MultShuff;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, log10(yval), '-', 'Color', plotcols{j});
    end
    axis tight
    ylim([-5 0]);
    line([predur_global predur_global], ylim);
    line(xlim, log10([alpha alpha]), 'Color', 'r', 'LineStyle', '--');
    
    % proportion branches with signficant p-val
    
    % -- 3) anova - r2
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.EtaSquared_shufflecollect;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, yval, '-', 'Color', plotcols{j});
    end
    axis tight
    ylim([-0.1 0.4])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
    
    
end



%% ====== COMBINE ALL THOSE CONTEXT METRICS ACROSS ALL BRANCHES FOR A GIVEN NEURON [ FIG = BRANCH]
alpha = 0.05; % for anova p-val

% figure out how many branches
% tmp = [];
% for nn=1:NumNeurons
%     tmp = [tmp ...
%         length(DATSTRUCT.neuron(nn).(branchtype).branch_master)];
% end
numBranches = length(allsyls);

plotcols = lt_make_plot_colors(NumNeurons, 0, 0);

for j=1:numBranches
    figcount=1;
    subplotrows=3;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    
    branchmaster = allsyls(j);
    
    %     plotcols = lt_make_plot_colors(numBranches, 0, 0);
    numNeuronsThisBranch = 0;
    % === ACTUAL DATA
    % -- 1) d-prime
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['dprime, BRANCH ' branchmaster]);
    ylabel('ACTUAL dat');
    xlabel('sec');
    for nn = 1:NumNeurons;
        if length(DATSTRUCT.neuron(nn).(branchtype).branch_master) < j
            continue
        end
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.dprime.xbins;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.dprime.meanpairwise;
        plot(xbins, yval, '-', 'Color', plotcols{nn});
        numNeuronsThisBranch = numNeuronsThisBranch +1;
    end
    axis tight
    ylim([-1 2])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
    
    % ----- IF ONLY ONE NEURON FOR THIS BRANCH, THEN CLOSE FIGURE AND
    % COTNINEU ETO NEXT BRANCH
    if numNeuronsThisBranch<2
        close(gcf)
        continue
    end
    
    % -- 2) anova pval
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('log10(p-val)');
    for nn = 1:NumNeurons;
        if length(DATSTRUCT.neuron(nn).(branchtype).branch_master) < j
            continue
        end
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.xbins;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.Pval;
        plot(xbins, log10(yval), '-', 'Color', plotcols{nn});
    end
    axis tight
    ylim([-5 0]);
    line([predur_global predur_global], ylim);
    line(xlim, log10([alpha alpha]), 'Color', 'r', 'LineStyle', '--');
    
    % proportion branches with signficant p-val
    
    % -- 3) anova - r2
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('eta^2');
    for nn = 1:NumNeurons;
        if length(DATSTRUCT.neuron(nn).(branchtype).branch_master) < j
            continue
        end
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.xbins;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.EtaSquared;
        plot(xbins, yval, '-', 'Color', plotcols{nn});
    end
    axis tight
    ylim([-0.1 0.4])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
    
    
    % ============= NEG CONTROL (CONTEXT SHUFFLED)
    % -- 4) d-prime
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('mean x-shuff');
    ylabel('NEG CTRL (CTXT SHUFFLE)');
    for nn = 1:NumNeurons;
        if length(DATSTRUCT.neuron(nn).(branchtype).branch_master) < j
            continue
        end
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.dprime.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.dprime.meanpairwise_MultShuff;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, yval, '-', 'Color', plotcols{nn});
    end
    axis tight
    ylim([-1 2])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
    
    
    % -- 2) anova pval
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('mean x-shuff');
    for nn = 1:NumNeurons;
        if length(DATSTRUCT.neuron(nn).(branchtype).branch_master) < j
            continue
        end
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.Pval_MultShuff;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, log10(yval), '-', 'Color', plotcols{nn});
    end
    axis tight
    ylim([-5 0]);
    line([predur_global predur_global], ylim);
    line(xlim, log10([alpha alpha]), 'Color', 'r', 'LineStyle', '--');
    
    % proportion branches with signficant p-val
    
    % -- 3) anova - r2
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    for nn = 1:NumNeurons;
        if length(DATSTRUCT.neuron(nn).(branchtype).branch_master) < j
            continue
        end
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.EtaSquared_shufflecollect;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, yval, '-', 'Color', plotcols{nn});
    end
    axis tight
    ylim([-0.1 0.4])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
    
    
    % ================= POS CONTROL (DIFF SYLS)
    % -- 4) d-prime
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('mean x-shuff');
    ylabel('POST CTRL (SLAVE SYls)');
    for nn = 1:NumNeurons;
        if length(DATSTRUCT.neuron(nn).(branchtype).branch_master) < j
            continue
        end
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.dprime.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.dprime.meanpairwise_MultShuff;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, yval, '-', 'Color', plotcols{nn});
    end
    axis tight
    ylim([-1 2])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
    
    
    % -- 2) anova pval
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('mean x-shuff');
    for nn = 1:NumNeurons;
        if length(DATSTRUCT.neuron(nn).(branchtype).branch_master) < j
            continue
        end
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.Pval_MultShuff;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, log10(yval), '-', 'Color', plotcols{nn});
    end
    axis tight
    ylim([-5 0]);
    line([predur_global predur_global], ylim);
    line(xlim, log10([alpha alpha]), 'Color', 'r', 'LineStyle', '--');
    
    % proportion branches with signficant p-val
    
    % -- 3) anova - r2
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    for nn = 1:NumNeurons;
        if length(DATSTRUCT.neuron(nn).(branchtype).branch_master) < j
            continue
        end
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.xbins_MultShuff;
        yval_cell = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.EtaSquared_shufflecollect;
        xmax = min(cellfun(@length, yval_cell));
        yval = [];
        for jj=1:length(yval_cell)
            yval = [yval; yval_cell{jj}(1:xmax)];
        end
        xbins = xbins(1:xmax);
        yval = mean(yval, 1);
        plot(xbins, yval, '-', 'Color', plotcols{nn});
    end
    axis tight
    ylim([-0.1 0.4])
    lt_plot_zeroline
    line([predur_global predur_global], ylim);
    
    
end

%% ====== COMBINE PLOT ACROSS ALL BRANCH POINTS AND ALL NEURONS [BRANCH X NEURON]

figcount=1;
subplotrows=3;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];

overlayRaw = 1;
overlaySummary = 0;
plotline = ':';

alpha = 0.05;

% ============= 1) DPRIME (actual)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('dprime');
ylabel('ACTUAL dat');
xlabel('sec');

Xbins_all={};
Yval_all={};

Presyl_onoff_all = [];
Postsyl_onoff_all = [];
Token_onoff_all= [];
for nn=1:NumNeurons
    
    numBranches = length(DATSTRUCT.neuron(nn).(branchtype).branch_master);
    
    % === ACTUAL DATA
    % -- 1) d-prime
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.dprime.xbins;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.dprime.meanpairwise;
        
        % ----- COLLECT
        Xbins_all = [Xbins_all xbins];
        Yval_all = [Yval_all yval];

        % --- COLLECT PRESYL TIMING
        DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.AllPresylOnOff;
        Presyl_onoff_all = [Presyl_onoff_all; ...
            DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.AllPresylOnOff];
        Postsyl_onoff_all = [Postsyl_onoff_all; ...
            DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.AllPostsylOnOff];
        Token_onoff_all = [Token_onoff_all; ...
            DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.AllTokenOnOffsets];
    end
end
% --- PLOT
xmax = min(cellfun(@length, Xbins_all));
xall = [];
yall = [];
for i=1:length(Yval_all)
    x=Xbins_all{i};
    y=Yval_all{i};
    
    if overlayRaw ==1
        plot(x, y, plotline, 'Color', [0.7 0.3 0.3]);
    end
    % -- collect all
    xall = [xall; x(1:xmax)];
    yall = [yall; y(1:xmax)];
end

% --- plot distributions at eacht ime point
if overlaySummary==1
    distributionPlot(yall, 'xValues', median(xall, 1), 'showMM', 0, 'color', [0.4 0.7 0.7])
end

axis tight
ylim([-1 2])
lt_plot_zeroline
line([predur_global predur_global], ylim);

% set(gca, 'XTick', 1:4:length(xall), 'XTickLabel', xtickvals);
% set(gca, 'XTick', 1:4:length(xall));
if strcmp(branchtype, 'conv')    
xmean = mean(Presyl_onoff_all, 1);
    line(xmean, [-0.1 -0.1], 'Color', [0.3 0.3 0.8], 'LineWidth', 2);
else
    xmean = mean(Postsyl_onoff_all, 1);
    line(xmean, [-0.1 -0.1], 'Color', [0.3 0.3 0.8], 'LineWidth', 2);
end
    xmean = mean(Token_onoff_all, 1);
    line(xmean, [-0.1 -0.1], 'Color', [0.3 0.3 0.8], 'LineWidth', 2);



% ============== 2) ANOVA PVAL
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('log10(pval)');
ylabel('ACTUAL dat');
xlabel('sec');

Xbins_all={};
Yval_all={};
for nn=1:NumNeurons
    
    numBranches = length(DATSTRUCT.neuron(nn).(branchtype).branch_master);
    
    % === ACTUAL DATA
    % -- 1) d-prime
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.xbins;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.Pval;
        
        % ----- COLLECT
        Xbins_all = [Xbins_all xbins];
        Yval_all = [Yval_all yval];
        
    end
end
% --- PLOT
xmax = min(cellfun(@length, Xbins_all));
xall = [];
yall = [];
for i=1:length(Yval_all)
    x=Xbins_all{i};
    y=Yval_all{i};
    if overlayRaw ==1
        
        plot(x, log10(y), plotline, 'Color', [0.7 0.3 0.3]);
    end
    % -- collect all
    xall = [xall; x(1:xmax)'];
    yall = [yall; y(1:xmax)];
end

% --- plot distributions at eacht ime point
if overlaySummary==1
    distributionPlot(log10(yall), 'xValues', median(xall, 1), 'showMM', 0, 'color', [0.4 0.7 0.7])
end
axis tight
ylim([-5 0]);
line([predur_global predur_global], ylim);
line(xlim, log10([alpha alpha]), 'Color', 'r', 'LineStyle', '--');
    
xmean = mean(Presyl_onoff_all, 1);
    line(xmean, [-4 -4], 'Color', [0.3 0.3 0.8], 'LineWidth', 2);


% ============ 2.5 - how many of the pvals (at each time point) are significant?
% - for a given alpha.
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('% pval<alpha');
ylabel('fraction');
xlabel('sec');

y_proportion = sum(yall<alpha, 1)./size(yall, 1);
lt_plot_bar(median(xall, 1), y_proportion);

axis tight;
line([predur_global predur_global], ylim);
ylim([0 0.6]);



% ============= 3) R2 anova (actual)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('eta^2');
ylabel('ACTUAL dat');
xlabel('sec');

Xbins_all={};
Yval_all={};
for nn=1:NumNeurons
    
    numBranches = length(DATSTRUCT.neuron(nn).(branchtype).branch_master);
    
    
    % === ACTUAL DATA
    % -- 1) d-prime
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.xbins;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat.anova.EtaSquared;
        
        % ----- COLLECT
        Xbins_all = [Xbins_all xbins];
        Yval_all = [Yval_all yval];
        
    end
end
% --- PLOT
xmax = min(cellfun(@length, Xbins_all));
xall = [];
yall = [];
for i=1:length(Yval_all)
    x=Xbins_all{i};
    y=Yval_all{i};
    
    if overlayRaw ==1
        plot(x, y, plotline, 'Color', [0.7 0.3 0.3]);
    end
    % -- collect all
    xall = [xall; x(1:xmax)'];
    yall = [yall; y(1:xmax)];
end

% --- plot distributions at eacht ime point
if overlaySummary==1
    distributionPlot(yall, 'xValues', median(xall, 1), 'showMM', 0, 'color', [0.4 0.7 0.7])
end
axis tight
ylim([-0.1 0.4])
lt_plot_zeroline
line([predur_global predur_global], ylim);

xmean = mean(Presyl_onoff_all, 1);
line(xmean, [-0.05 -0.05], 'Color', [0.3 0.3 0.8], 'LineWidth', 2);




% ============= 4) d - prime (NEG CONTROL)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean x-shuff');
ylabel('NEG CTRL (CTXT SHUFFLE)');

Xbins_all={};
Yval_all={};
for nn=1:NumNeurons
    
    numBranches = length(DATSTRUCT.neuron(nn).(branchtype).branch_master);
    
    
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.dprime.xbins_MultShuff;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.dprime.meanpairwise_MultShuff;
        
        % take mean yval
        xtmp = min(cellfun(@length, yval));
        tmptmp = [];
        for jj=1:length(yval)
            tmptmp = [tmptmp; yval{jj}(1:xtmp)];
        end
        yval = mean(tmptmp, 1);
        xbins = xbins(1:xtmp);
        
        % ----- COLLECT
        Xbins_all = [Xbins_all xbins];
        Yval_all = [Yval_all yval];
        
    end
end
% --- PLOT
xmax = min(cellfun(@length, Yval_all));
xall = [];
yall = [];
for i=1:length(Yval_all)
    y=Yval_all{i};
    x=Xbins_all{i};
    
    if overlayRaw ==1
        plot(x, y, plotline, 'Color', [0.7 0.3 0.3]);
    end
    % -- collect all
    xall = [xall; x(1:xmax)];
    yall = [yall; y(1:xmax)];
end

% --- plot distributions at eacht ime point
if overlaySummary==1
    distributionPlot(yall, 'xValues', median(xall, 1), 'showMM', 0, 'color', [0.4 0.7 0.7])
end
axis tight
ylim([-1 2])
lt_plot_zeroline
line([predur_global predur_global], ylim);



% ============= 5) pval (anova) (NEG CONTROL)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean x-shuff');

Xbins_all={};
Yval_all={};
for nn=1:NumNeurons
    
    numBranches = length(DATSTRUCT.neuron(nn).(branchtype).branch_master);
    
    
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.xbins_MultShuff;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.Pval_MultShuff;
        
        % take mean yval
        xtmp = min(cellfun(@length, yval));
        tmptmp = [];
        for jj=1:length(yval)
            tmptmp = [tmptmp; yval{jj}(1:xtmp)];
        end
        yval = mean(tmptmp, 1);
        xbins = xbins(1:xtmp);
        
        % ----- COLLECT
        Xbins_all = [Xbins_all xbins];
        Yval_all = [Yval_all yval];
        
    end
end
% --- PLOT
xmax = min(cellfun(@length, Yval_all));
xall = [];
yall = [];
for i=1:length(Yval_all)
    y=Yval_all{i};
    x=Xbins_all{i};
    if overlayRaw ==1
        plot(x, log10(y), plotline, 'Color', [0.7 0.3 0.3]);
    end
    
    % -- collect all
    xall = [xall; x(1:xmax)'];
    yall = [yall; y(1:xmax)];
end

% --- plot distributions at eacht ime point
if overlaySummary==1
    distributionPlot(log10(yall), 'xValues', median(xall, 1), 'showMM', 0, 'color', [0.4 0.7 0.7])
end
axis tight
ylim([-5 0]);
line([predur_global predur_global], ylim);
line(xlim, log10([alpha alpha]), 'Color', 'r', 'LineStyle', '--');



% ============ 5.5 - how many of the pvals (at each time point) are significant?
% - for a given alpha.
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('% pval<alpha');
ylabel('fraction');
xlabel('sec');

y_proportion = sum(yall<alpha, 1)./size(yall, 1);
lt_plot_bar(median(xall, 1), y_proportion);

axis tight;
line([predur_global predur_global], ylim);
ylim([0 0.6]);



% ============= 6) Eta^2 (anova) (NEG CONTROL)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

Xbins_all={};
Yval_all={};
for nn=1:NumNeurons
    
    numBranches = length(DATSTRUCT.neuron(nn).(branchtype).branch_master);
    
    
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.xbins_MultShuff;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).Control_ContextShuff.anova.EtaSquared_shufflecollect;
        
        % take mean yval
        xtmp = min(cellfun(@length, yval));
        tmptmp = [];
        for jj=1:length(yval)
            tmptmp = [tmptmp; yval{jj}(1:xtmp)];
        end
        yval = mean(tmptmp, 1);
        xbins = xbins(1:xtmp);
        
        % ----- COLLECT
        Xbins_all = [Xbins_all xbins];
        Yval_all = [Yval_all yval];
        
    end
end
% --- PLOT
xmax = min(cellfun(@length, Yval_all));
xall = [];
yall = [];
for i=1:length(Yval_all)
    y=Yval_all{i};
    x=Xbins_all{i};
    
    if overlayRaw ==1
        plot(x, y, plotline, 'Color', [0.7 0.3 0.3]);
    end
    
    % -- collect all
    xall = [xall; x(1:xmax)'];
    yall = [yall; y(1:xmax)];
end

% --- plot distributions at eacht ime point
if overlaySummary==1
    distributionPlot(yall, 'xValues', median(xall, 1), 'showMM', 0, 'color', [0.4 0.7 0.7])
end
axis tight
ylim([-0.1 0.4])
lt_plot_zeroline
line([predur_global predur_global], ylim);




% ============= 7) d - prime (POS CONTROL)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean x-shuff');
ylabel('POST CTRL (SLAVE SYls)');

Xbins_all={};
Yval_all={};
for nn=1:NumNeurons
    
    numBranches = length(DATSTRUCT.neuron(nn).(branchtype).branch_master);
    
    
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.dprime.xbins_MultShuff;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.dprime.meanpairwise_MultShuff;
        
        % take mean yval
        xtmp = min(cellfun(@length, yval));
        tmptmp = [];
        for jj=1:length(yval)
            tmptmp = [tmptmp; yval{jj}(1:xtmp)];
        end
        yval = mean(tmptmp, 1);
        xbins = xbins(1:xtmp);
        
        % ----- COLLECT
        Xbins_all = [Xbins_all xbins];
        Yval_all = [Yval_all yval];
        
    end
end
% --- PLOT
xmax = min(cellfun(@length, Yval_all));
xall = [];
yall = [];
for i=1:length(Yval_all)
    y=Yval_all{i};
    x=Xbins_all{i};
    
    if overlayRaw ==1
        plot(x, y, plotline, 'Color', [0.7 0.3 0.3]);
    end
    
    % -- collect all
    xall = [xall; x(1:xmax)];
    yall = [yall; y(1:xmax)];
end

% --- plot distributions at eacht ime point
if overlaySummary==1
    distributionPlot(yall, 'xValues', median(xall, 1), 'showMM', 0, 'color', [0.4 0.7 0.7])
end

axis tight
ylim([-1 2])
lt_plot_zeroline
line([predur_global predur_global], ylim);


% ============= 8) d - prime (POS CONTROL)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean x-shuff');

Xbins_all={};
Yval_all={};
for nn=1:NumNeurons
    
    numBranches = length(DATSTRUCT.neuron(nn).(branchtype).branch_master);
    
    
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.xbins_MultShuff;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.Pval_MultShuff;
        
        % take mean yval
        xtmp = min(cellfun(@length, yval));
        tmptmp = [];
        for jj=1:length(yval)
            tmptmp = [tmptmp; yval{jj}(1:xtmp)];
        end
        yval = mean(tmptmp, 1);
        xbins = xbins(1:xtmp);
        
        % ----- COLLECT
        Xbins_all = [Xbins_all xbins];
        Yval_all = [Yval_all yval];
        
    end
end
% --- PLOT
xmax = min(cellfun(@length, Yval_all));
xall = [];
yall = [];
for i=1:length(Yval_all)
    y=Yval_all{i};
    x=Xbins_all{i};
    
    if overlayRaw ==1
        plot(x, log10(y), plotline, 'Color', [0.7 0.3 0.3]);
    end
    
    % -- collect all
    xall = [xall; x(1:xmax)'];
    yall = [yall; y(1:xmax)];
end

% --- plot distributions at eacht ime point
if overlaySummary==1
    distributionPlot(log10(yall), 'xValues', median(xall, 1), 'showMM', 0, 'color', [0.4 0.7 0.7])
end
axis tight
ylim([-5 0]);
line([predur_global predur_global], ylim);
line(xlim, log10([alpha alpha]), 'Color', 'r', 'LineStyle', '--');



% ============ 8.5 - how many of the pvals (at each time point) are significant?
% - for a given alpha.
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('% pval<alpha');
ylabel('fraction');
xlabel('sec');

y_proportion = sum(yall<alpha, 1)./size(yall, 1);
lt_plot_bar(median(xall, 1), y_proportion);

axis tight;
line([predur_global predur_global], ylim);
ylim([0 0.6]);



% ============= 9) anova - r2 (POS CONTROL)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

Xbins_all={};
Yval_all={};
for nn=1:NumNeurons
    
    numBranches = length(DATSTRUCT.neuron(nn).(branchtype).branch_master);
    
    
    for j=1:numBranches
        if isempty(DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ActualDat)
            continue
        end
        
        xbins = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.xbins_MultShuff;
        yval = DATSTRUCT.neuron(nn).(branchtype).branch_master(j).ControlPos_SlaveSyls.anova.EtaSquared_shufflecollect;
        
        % take mean yval
        xtmp = min(cellfun(@length, yval));
        tmptmp = [];
        for jj=1:length(yval)
            tmptmp = [tmptmp; yval{jj}(1:xtmp)];
        end
        yval = mean(tmptmp, 1);
        xbins = xbins(1:xtmp);
        
        % ----- COLLECT
        Xbins_all = [Xbins_all xbins];
        Yval_all = [Yval_all yval];
        
    end
end
% --- PLOT
xmax = min(cellfun(@length, Yval_all));
xall = [];
yall = [];
for i=1:length(Yval_all)
    y=Yval_all{i};
    x=Xbins_all{i};
    
    if overlayRaw ==1
        plot(x, y, plotline, 'Color', [0.7 0.3 0.3]);
    end
    
    % -- collect all
    xall = [xall; x(1:xmax)'];
    yall = [yall; y(1:xmax)];
end

% --- plot distributions at eacht ime point
if overlaySummary==1
    distributionPlot(yall, 'xValues', median(xall, 1), 'showMM', 0, 'color', [0.4 0.7 0.7])
end
axis tight
ylim([-0.1 0.4])
lt_plot_zeroline
line([predur_global predur_global], ylim);




%% OBSOLETE - REPLACED MY METHOD ABOVE, IN WHICH I FIERST COLLECT DATA, THEN
% PLOT....
%% ===== for each branch point, plot every neuron's response [CONV]
% if (0)
%
% hallfigs=findobj;
% maxfig = max(hallfigs);
%
% for nn=1:NumNeurons
%
%     % 1) -- load dat
%     cd(NeuronDatabase.global.basedir);
%
%     % - find day folder
%     dirdate=NeuronDatabase.neurons(nn).date;
%     tmp=dir([dirdate '*']);
%     assert(length(tmp)==1, 'PROBLEM - issue finding day folder');
%     cd(tmp(1).name);
%
%     % - load data for this neuron
%     batchf=NeuronDatabase.neurons(nn).batchfile;
%     channel_board=NeuronDatabase.neurons(nn).chan;
%     [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);
%
%     hsplots = cell(length(allsyls), 1);
%     for i=1:length(allsyls)
%         syl_master = allsyls(i); % this syl defines the brainch point - i..e if conv, then is syl2
%
%         figure(i+maxfig); hold on;
%         lt_plot_format;
%         hsplot = lt_subplot(8, 3, nn); hold on;
%         title(['neuron ' num2str(nn)]);
%         hsplots{i} = [hsplots{i} hsplot];
%
%         plotcols = lt_make_plot_colors(length(allsyls), 0, 0); % prepare for 6 transitions
%
%         % === for this neuron and branch point, get all transitions and
%         % plot each one.
%         for ii=1:length(allsyls)
%             syl_slave = allsyls(ii);
%
%             if strcmp(branchtype, 'conv')
%                 regexpr_str = [syl_slave '(' syl_master ')'];
%                 alignByOnset = 1;
%                 predur = 0.15;
%                 postdur = 0;
%             elseif strcmp(branchtype, 'div')
%                 regexpr_str = ['(' syl_master ')' syl_slave];
%                 alignByOnset = 0;
%                 predur = 0.1;
%                 postdur = 0.05;
%             end
%
%             plotcol=plotcols{ii}; % color for this transition
%
%             % ----- criteria to skip
%             if strcmp(syl_slave, '-') || strcmp(syl_master, '-');
%                 continue
%             end
%             % ------
%
%             % ================ GET SMOOTHED FIRING RATE FOR THIS MOTIF
%             % 2) -- extract dat
%             [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
%                 regexpr_str, predur, postdur, alignByOnset, '');
%
%             % -- throw out if not enough datapoints
%             if length(SegmentsExtract) < minNperTran
%                 continue
%             end
%
%             % 3) ---------------------- get smoothed firing rate over all trials
%             numtrials=length(SegmentsExtract);
%             clustnum=NeuronDatabase.neurons(nn).clustnum;
%             Yspks={};
%             for k=1:numtrials
%                 inds=SegmentsExtract(k).spk_Clust==clustnum;
%                 spktimes=SegmentsExtract(k).(spktimefield)(inds);
%                 Yspks{k}=spktimes;
%             end
%
%             % -- convert to smoothed rate
%             [xbin, ~, ~, ymean_hz, ysem_hz] = lt_neural_plotRastMean(Yspks, window, windshift, 0, '');
%
%
%             % ===================== PLOT FOR THIS NEURON
%             shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcol}, 1);
%             lt_plot_text(xbin(end)+0.01, ymean_hz(end), syl_slave, plotcol);
%
%             % ------- final things
%             line([predur predur], ylim, 'Color', 'k');
%             axis tight
%             lt_plot_zeroline;
%             tmp=get(gca, 'Ylim');
%             set(gca, 'Ylim', [-16 tmp(2)]);
%
%             % --- put syl duration lines.
%             if strcmp(branchtype, 'conv')
%                 token_offtime = median(SongDat.AllOffsets([SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]) ...
%                     - [SegmentsExtract.global_ontime_motifInclFlank]);
%                 line([predur token_offtime], -rand(1)*10 - [5 5], 'Color', plotcol, 'LineWidth', 1.5)
%
%             elseif strcmp(branchtype, 'div')
%                 token_ontime = median(SongDat.AllOnsets([SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]) ...
%                     - [SegmentsExtract.global_ontime_motifInclFlank]);
%                 line([token_ontime predur], -rand(1)*10 - [5 5], 'Color', plotcol, 'LineWidth', 1.5)
%             end
%
%
%         end
%
%         linkaxes(hsplots{i}, 'x');
%         if nn == NumNeurons
%             if strcmp(branchtype, 'conv')
%                 lt_subtitle(['CONV to ' syl_master]);
%             elseif strcmp(branchtype, 'div');
%                 lt_subtitle(['DIV from ' syl_master]);
%             end
%         end
%     end
% end
%
% end
%


end

function         [FRdiffAcrossBranch, XbinsAcrossBranch, FRdiff_thisbranch, Xbin_forplot_thisbranch]  ...
    = fn_get_diff_from_slave_contour(AllSlaveContours, plot_on, AllXbins, ...
    AllTokenOnOffsets, FRdiffAcrossBranch, XbinsAcrossBranch, diffmetric, AllSlaveContours_std)

if ~exist('diffmetric', 'var');
    diffmetric = 'frdiff'; % default is frate diff.
end

xmax = min(cellfun(@length, AllSlaveContours));
% 1) get all paired differences of contour means
numcontours = length(AllSlaveContours);
if numcontours > 1
    tmp_alldiffs = [];
    for k=1:numcontours-1
        for kk=k+1:numcontours
            if strcmp(diffmetric, 'frdiff')
                tmp_diff = abs(AllSlaveContours{k}(1:xmax) - AllSlaveContours{kk}(1:xmax));
            elseif strcmp(diffmetric, 'dprime') % takes absolute d-prime (absolute because which contorur subtract which is arbitrary)
                numerator = AllSlaveContours{k}(1:xmax) - AllSlaveContours{kk}(1:xmax);
                denominator = sqrt(0.5 * (AllSlaveContours_std{k}(1:xmax).^2 + ...
                    AllSlaveContours_std{kk}(1:xmax).^2));
                dprime = numerator./denominator;
                tmp_diff = abs(dprime);
            end
            
            tmp_alldiffs = [tmp_alldiffs; tmp_diff];
        end
    end
    
    if plot_on ==1
        if strcmp(diffmetric, 'dprime')
            multiplier = 1;
        else
            multiplier = 1;
        end
        % ====== PLOT FOR THIS NEURON
        plot(AllXbins{1}(1:xmax), mean(tmp_alldiffs, 1) * multiplier, '-ok'); % this contains all diffs for this branch
        
        % === PLOT LINE FOR SYL
        if strcmp(diffmetric, 'dprime');
            line(mean(AllTokenOnOffsets,1), -rand(1)*0.3 - [0.6 0.6], 'Color', [0.2 0.7 0.7], 'LineWidth', 1.5)
            ylim([-1 2]);
        else
            line(mean(AllTokenOnOffsets,1), -rand(1)*10 - [5 5], 'Color', [0.2 0.7 0.7], 'LineWidth', 1.5)
        end
        line([mean(AllTokenOnOffsets(:,1)) mean(AllTokenOnOffsets(:,1))], ylim, 'Color', 'k');
        
        axis tight;
        lt_plot_zeroline;
        
        if strcmp(diffmetric, 'dprime')
            ylabel(['dprime (x' num2str(multiplier) ')']);
        elseif strcmp(diffmetric, 'frdiff');
            ylabel('fr diff(hz)');
        end
        
    end
    
    % ==== collect
    FRdiff_thisbranch = mean(tmp_alldiffs,1);
    Xbin_forplot_thisbranch = AllXbins{1}(1:xmax);
    FRdiffAcrossBranch = [FRdiffAcrossBranch mean(tmp_alldiffs, 1)];
    XbinsAcrossBranch = [XbinsAcrossBranch AllXbins{1}];
    
end


end

function    fn_plot_mean_diff_for_neuron(XbinsAcrossBranch, FRdiffAcrossBranch, plotcol)


x2 = max(cellfun(@length, FRdiffAcrossBranch)); % not controlling for sample size
tmptmp=nan(length(FRdiffAcrossBranch), x2);
for lll = 1:length(FRdiffAcrossBranch)
    tmptmp(lll, 1:length(FRdiffAcrossBranch{lll})) = FRdiffAcrossBranch{lll};
end
xdiff = XbinsAcrossBranch{1}(2)-XbinsAcrossBranch{1}(1);
xbins=XbinsAcrossBranch{1}(1):xdiff:(x2-1)*xdiff+XbinsAcrossBranch{1}(1);
%     plot(xbins, nanmean(tmptmp,1), '-k');
if length(FRdiffAcrossBranch)>1
    shadedErrorBar(xbins, nanmean(tmptmp, 1), lt_sem(tmptmp), { 'Color', plotcol},1 );
else
    plot(xbins, tmptmp, '-s', 'Color', plotcol);
end

end