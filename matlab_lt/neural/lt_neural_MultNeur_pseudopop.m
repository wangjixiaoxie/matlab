function lt_neural_MultNeur_pseudopop(NeuronDatabase, TRANSMATRIX, branchtype)


minNeuronsPerBranch = 1; % only plot a branch point if N neurons have this labeled
minNperTran = 8; % min sample size
motif_predur=0.13;
motif_postdur=0.05;

spktimefield='spk_Times';
window = 0.02;
windshift = 0.004;

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
                postdur = 0.05;
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
                %
            elseif strcmp(branchtype, 'div')
                token_ontime = median(SongDat.AllOnsets([SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]) ...
                    - [SegmentsExtract.global_ontime_motifInclFlank]);
                %                 line([token_ontime predur], -rand(1)*10 - [5 5], 'Color', plotcol, 'LineWidth', 1.5)
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


%% ========== PLOT FIRIGN RATES AT ALL BRANCHES

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
            
            if ~isempty(xbin)
                shadedErrorBar(xbin, ymean_hz, ysem_hz, {'Color', plotcols{ii}}, 1);
                lt_plot_text(xbin(end)+0.01, ymean_hz(end), syl_slave, plotcols{ii});
                line(token_onset_offset, -rand(1)*10 - [5 5], 'Color', plotcols{ii}, 'LineWidth', 1.5)
                %                     end
            end
        end
        
        axis tight
        
        % ---- Overlay vert line
        predur = DATSTRUCT.neuron(nn).(branchtype).params.predur;
        line([predur predur], ylim, 'Color', 'k');
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

%% ======== FOR EACH BRANCH POINT, GET RUNNING AVERAGE OF FIRING RATE DEVIATIONS FOR ALL PAIRS OF TRANSITIONS
% =================== ACTUAL, GET ALL PAIRWISE D-PRIMES (BETWEEN ALL
% TRANSITIONS FOR A CERTAIN BRANCH POINT

for i=1:length(allsyls)
    syl_master = allsyls(i);
    if strcmp(syl_master, '-')
        continue
    end
    lt_figure; hold on;
    
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
        xmax = min(cellfun(@length, AllSlaveContours));
        % 1) get all paired differences of contour means
        numcontours = length(AllSlaveContours);
        if numcontours > 1
            tmp_alldiffs = [];
            for k=1:numcontours-1
                for kk=k+1:numcontours
                    tmp_diff = abs(AllSlaveContours{k}(1:xmax) - AllSlaveContours{kk}(1:xmax));
                    tmp_alldiffs = [tmp_alldiffs; tmp_diff];
                end
            end
            
            % - PLOT
            plot(AllXbins{1}(1:xmax), mean(tmp_alldiffs, 1), '-o');
            
            % ---- OVERLAY SYL ONSET OFFSET
            line(mean(AllTokenOnOffsets,1), -rand(1)*10 - [5 5], 'Color', 'k', 'LineWidth', 1.5)
            axis tight
            
            % ---- Overlay vert line
            predur = DATSTRUCT.neuron(nn).(branchtype).params.predur;
            line([predur predur], ylim, 'Color', 'k');
            
            lt_plot_zeroline;
            
        end
    end
    
    if strcmp(branchtype, 'conv')
        lt_subtitle(['CONV to ' syl_master]);
    elseif strcmp(branchtype, 'div');
        lt_subtitle(['DIV from ' syl_master]);
    end
end


%% ================= PLOT ALL DIFFS FOR EACH NEURON (I.E. EACH NEURON, ONE PLOT)
% COLLECT AND PLOT FOR EACH BRANCH POINT

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

%% ========== TROUBLESHOOTING ABOVE - 
% for a given neuron, and given master syl, plot 1) all actual transitions,
% 2) negative control - i.e. all shuffled, 3) positive control - i.e. all
% diff branches

syl_master_wanted = 'h';
hsplots = [];


for nn=1:NumNeurons
figcount=1;
subplotrows=3;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];

    ylabel('mean FR diff (hz)');
    xlabel('ms');
    
    
    % ================================ 1) ALL ACTUAL TRANSITIONS
    FRdiffAcrossBranch = {};
    XbinsAcrossBranch = {};
    for i= strfind(allsyls, syl_master_wanted);
%     for i= [2 6];
        syl_master = allsyls(i);
        if strcmp(syl_master, '-')
            continue
        end
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);       
        title(['syl master: ' syl_master]);
        hsplots = [hsplots hsplot];
        
        % for this master (i.e. this branch) collect all contours, then
        % take difference between those contorus.
        AllSlaveContours = {};
        AllSlaveContours_std = {};
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
            
            % --- PLOT
            plotcol = [rand rand rand];
            plot(xbin, ymean_hz, '-', 'Color', plotcol);
            lt_plot_text(xbin(end), ymean_hz(end), allsyls(ii), plotcol);
            line(token_onset_offset,-5+2*rand, 'Color', plotcol)
            
            % --- SAVE
            AllXbins = [AllXbins xbin];
            AllSlaveContours = [AllSlaveContours ymean_hz];
            AllSlaveContours_std = [AllSlaveContours_std ystd_hz];
            AllTokenOnOffsets = [AllTokenOnOffsets; token_onset_offset];
            %             end
        end
        
        % ============= DO STUFF WITH THE SLAVE CONTOURS
        plot_on = 1; % plots for a given neuron, mean for this diff metric
        diffmetric = 'dprime';
%         diffmetric = 'frdiff';
        [FRdiffAcrossBranch, XbinsAcrossBranch] = fn_get_diff_from_slave_contour(...
            AllSlaveContours, plot_on, AllXbins, AllTokenOnOffsets, FRdiffAcrossBranch, ...
            XbinsAcrossBranch, diffmetric, AllSlaveContours_std);
        
        
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
    if length(FRdiffAcrossBranch)>1
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
    end
    
    
    % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++===
    % ++++++++++++++++++++++++++++++++++++++++++++++ CONTROLS
    
    % ============ 1) for this branch, shuffle transitions - i.e.
    % assuming there is no difference in firing across transitions,
    % then what types of contours would we expect?
    FRdiffAcrossBranch = {};
    XbinsAcrossBranch = {};
    for i = strfind(allsyls, syl_master_wanted);
        disp(['master syl ' num2str(i)]);
        syl_master = allsyls(i);
        if strcmp(syl_master, '-')
            continue
        end
        
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);       
        title(['[shuff transitions: syl master: ' syl_master]);
hsplots = [hsplots hsplot];

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
            
            % == plot
            plotcol = [rand rand rand];
            plot(xbin, ymean_hz, '-', 'Color', plotcol);
            
            
            % -- save mean firing rate
            AllXbins = [AllXbins xbin];
            AllSlaveContours = [AllSlaveContours ymean_hz];
        end
        
        % ============= DO STUFF WITH THE SLAVE CONTOURS
        plot_on = 1; % plots for a given neuron, mean for this diff metric
        diffmetric = 'dprime';
%         diffmetric = 'frdiff';
        [FRdiffAcrossBranch, XbinsAcrossBranch] = fn_get_diff_from_slave_contour(...
            AllSlaveContours, plot_on, AllXbins, AllTokenOnOffsets, FRdiffAcrossBranch, ...
            XbinsAcrossBranch, diffmetric, AllSlaveContours_std);
        
    end
    
    % ---- PLOT MEAN DIFF (across branches) for this neuron
    if length(FRdiffAcrossBranch) >1
    plotcol = [0.2 0.9 0.2];
    fn_plot_mean_diff_for_neuron(XbinsAcrossBranch, FRdiffAcrossBranch, plotcol)
    % ====================================================================
    end
    
    lt_plot_zeroline
    
    
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
    for i = strfind(allsyls, syl_master_wanted);

        syl_master = allsyls(i);
        if strcmp(syl_master, '-')
            continue
        end
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);       
        title('diff syls control');
hsplots = [hsplots hsplot];

        
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
            
            % === plot
            plotcol = [rand rand rand];
            plot(xbin, ymean_hz, '-', 'Color', plotcol);
            lt_plot_text(xbin(end), ymean_hz(end), allsyls(ii), plotcol);
            
            % -- output
            AllSlaveContours = [AllSlaveContours ymean_hz];
            AllXbins = [AllXbins xbin];
        end
        
        % ---- Do stuff with slave contours
        plot_on = 1;
        diffmetric = 'dprime';
%         diffmetric = 'frdiff';
        [FRdiffAcrossBranch, XbinsAcrossBranch] = fn_get_diff_from_slave_contour(...
            AllSlaveContours, plot_on, AllXbins, AllTokenOnOffsets, FRdiffAcrossBranch, ...
            XbinsAcrossBranch, diffmetric, AllSlaveContours_std);
    end
    
    
    % ---- PLOT MEAN DIFF (across branches) for this neuron
    if length(FRdiffAcrossBranch) > 1
    plotcol = 'r';
    fn_plot_mean_diff_for_neuron(XbinsAcrossBranch, FRdiffAcrossBranch, plotcol);
    % ====================================================================
    end
    
    lt_plot_zeroline;
end

linkaxes(hsplots, 'xy');

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

function         [FRdiffAcrossBranch, XbinsAcrossBranch] = fn_get_diff_from_slave_contour(...
    AllSlaveContours, plot_on, AllXbins, AllTokenOnOffsets, FRdiffAcrossBranch, ...
    XbinsAcrossBranch, diffmetric, AllSlaveContours_std)

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
            multiplier = 50;
        else
            multiplier = 1;
        end
        % ====== PLOT FOR THIS NEURON
        plot(AllXbins{1}(1:xmax), mean(tmp_alldiffs, 1) * multiplier, '-o', 'Color', [0.2 0.7 0.7]); % this contains all diffs for this branch
        
        % === PLOT LINE FOR SYL
        line(mean(AllTokenOnOffsets,1), -rand(1)*10 - [5 5], 'Color', [0.2 0.7 0.7], 'LineWidth', 1.5)
    end
    
    % ==== collect
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