function lt_neural_MultNeur_MotifRasters_LearnSum(NeuronDatabase, motif_regexpr_str, motif_predur, ...
    motif_postdur, LinScaleGlobal, FFparams, OnlyPlotNoHit, UseEntireBaseline, TypesToPlot,TrialBinSize )

%% for each motif and each neuron, make one plot, showing progression over time of
% motif_regexpr_str={'g(h)'};
% motif_predur=0.2;
% motif_postdur=0.1;
% LinScaleGlobal=0; % 0:NONE; 1: global (across neurosn and motifs); 2: local (specific to neuron x motif)
%
%
% FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will
% FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
%     % +1 is 1 after token
% FFparams.FF_sylName='h'; % Optional: what syl do you expect this to be? if incompatible will raise error
%     % not required (can set as []);
% FFparams.cell_of_freqwinds={'h', [1100 2600], 'b', [2400 3500], ...
%             'v', [2450 4300]};
% % FFparams.cell_of_FFtimebins={'h', [0.042 0.058], 'b', [0.053 0.07], ...
% %             'v', [0.052 0.07]}; % in sec, relative to onset (i.e. see vector T)
% FFparams.cell_of_FFtimebins={'h', [0.034 0.038], 'b', [0.053 0.07], ...
%             'v', [0.052 0.07]}; % WN on g H
% % NOTE: will also determine whether was hit or miss, based on WN sound
% % detection.
%
% % LEARNING PARAMS
% WNchangeDateStrings={'05Oct2016-1348'};
%
% OnlyPlotNoHit=1; % then only plots trials that were not hit (WN)
% TypesToPlot={'mean', 'abs', 'std'}; % for d-prime summary metric


%% PARAMS
NumMotifs=length(motif_regexpr_str);
NumNeurons=length(NeuronDatabase.neurons);

% smoothing window (neural)
window=0.02; windshift=0.004;

plotTimesOnRaster=1; % hh:mm on raster plot

determineTrialBinBasedOnBaselineN=1; % if 1, then chooses binsize based on
% baseline num songs (i.e. baseNumSongs/DivisorBaseSongs); if 0, then = TrialBinSize
% note: assumes WNchangeDateStrings{1} is transition from baseline to WN
% on.
DivisorBaseSongs=2;
% TrialBinSize=10;

if LinScaleGlobal ==0
    spktimefield='spk_Times';
end

% for taking d-prime, relative to token onset
premotor_wind=[-0.08 0.02]; % in sec, relative to onset of token in motif.

% UseEntireBaseline = 0; % if 1, uses entire baseline, otherwise uses 1st bin.



%% EXTRACT DATA (for each neuron x motif)
MOTIFSTATS=struct;
for i=1:NumNeurons
    cd(NeuronDatabase.global.basedir);
    
    % - find day folder
    dirdate=NeuronDatabase.neurons(i).date;
    tmp=dir([dirdate '*']);
    assert(length(tmp)==1, 'PROBLEM - issue finding day folder');
    cd(tmp(1).name);
    
    % - load data for this neuron
    batchf=NeuronDatabase.neurons(i).batchfile;
    channel_board=NeuronDatabase.neurons(i).chan;
    [SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);
    
    % --- EXTRACT DAT
    % - do this one time for each desired motif
    for j=1:NumMotifs
        regexpr_str=motif_regexpr_str{j};
        predur=motif_predur; % sec
        postdur=motif_postdur; % sec
        alignByOnset=1;
        WHOLEBOUTS_edgedur=''; % OPTIONAL (only works if regexpr_str='WHOLEBOUTS', only keeps
        % those motifs that have long enough pre and post - LEAVE EMPTY TO GET ALL BOUTS
        [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
            regexpr_str, predur, postdur, alignByOnset, WHOLEBOUTS_edgedur, FFparams);
        
        
        % -------- SAVE TIMING OF SPIKES FOR THIS NEURON
        MOTIFSTATS.neurons(i).motif(j).SegmentsExtract=SegmentsExtract;
        MOTIFSTATS.neurons(i).motif(j).Params=Params;
    end
end

% pause;
% disp('PRESS ANYTHING TO CONTINUE');
% close all;


%% === MAKE SURE ALL TRIALS ARE IN TEMPORAL ORDER

for i=1:NumNeurons
    for m=1:NumMotifs
        
        segextract = MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
        
        all_datenums=[segextract.song_datenum];
        
        [~, inds] = sort(all_datenums);
        
        if all(diff(inds))~=1
            disp('HAD TO REORDER !! ---- not a problem');
            
            MOTIFSTATS.neurons(i).motif(m).SegmentsExtract=segextract(inds);
        end
    end
end


%% === SUMMARIZE BY TAKING MEAN DPRIME OVER PREMOTOR WINODW


for kkk=1:length(TypesToPlot)
    typetoplot=TypesToPlot{kkk};
    
    for i=1:NumNeurons
        figcount=1;
        subplotrows=8;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots1=[];
        hsplots2=[];
        
        clear WNchangeDateStrings;
WNchangeDateStrings{1}=NeuronDatabase.neurons(i).LEARN_WNonDatestr;
WNchangeDateStrings=[WNchangeDateStrings NeuronDatabase.neurons(i).LEARN_WNotherImportantDates];

for m=1:NumMotifs
            
            
            segextract=MOTIFSTATS.neurons(i).motif(m).SegmentsExtract;
            clustnum=NeuronDatabase.neurons(i).clustnum;
            
            
            % --- IF ONLY PLOT NO HIT, THEN FIRST PULL OUT THOSE TRIALS
            if OnlyPlotNoHit==1
                indsWithNoHits=[segextract.hit_WN]~=1;
                segextract=segextract(indsWithNoHits);
            end
            
            
            % ===================== collect spikes for all trials
            numtrials=length(segextract);
            Yspks={};
            for j=1:numtrials
                inds=segextract(j).spk_Clust==clustnum;
                spktimes=segextract(j).(spktimefield)(inds);
                Yspks{j}=spktimes;
            end
            
            % ===================== BASELINE mean and std spiking
            ymean_hz_base=[]; ystd_hz_base=[];
            if UseEntireBaseline==1
                datestring=WNchangeDateStrings{1}; % assumes this is WN on.
                dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                inds=[segextract.song_datenum] < dnum; %
                [~, ~, ~, ymean_hz_base, ~, ~, ystd_hz_base] = ...
                    lt_neural_plotRastMean(Yspks(inds), window, windshift, 0, '');
            elseif UseEntireBaseline==0
                % then will just use first bin as "baseline"
                inds=1:TrialBinSize;
                [~, ~, ~, ymean_hz_base, ~, ~, ystd_hz_base] = ...
                    lt_neural_plotRastMean(Yspks(inds), window, windshift, 0, '');
            end
            
            % ========================= GET RUNNING D-PRIME OVER TIME, IN
            % SPECIFIC PREMOTOR WINDOW
            DprimeVals=[];
            AbsDprimeVals=[];
            DprimeStd=[]; % variance of d-prime values across time bins
            CorrCoeff = []; 
            MedianDatenumVals = [];
            MedianDayVals=[];
            FFVals=[];
            AreAllTrialsBaseline = [];
            AreAllTrialsDurWN = [];
            
            for j=1:(numtrials-TrialBinSize+1);
                trialbins=j:(j+TrialBinSize-1);
                
                [x, ~, ~, ymean_hz, ~, ~, ystd_hz] =  ...
                    lt_neural_plotRastMean(Yspks(trialbins), window, windshift, 0, '');
                
                % -- calc d-prime
                xinds=1:min([length(ymean_hz) length(ymean_hz_base)]);
                numer_tmp=ymean_hz(xinds)-ymean_hz_base(xinds);
                denom_tmp=sqrt(0.5*(ystd_hz_base(xinds).^2 + ystd_hz(xinds).^2));
                dprime=numer_tmp./denom_tmp;
                
                % -- slice dprime within premotor window
                xinds_to_keep = x>=(motif_predur + premotor_wind(1) + window/2) ...
                    & x<=(motif_predur + premotor_wind(2) - window/2); % accounting for smoothing window size.
                
                dprime_val = mean(dprime(xinds_to_keep));
                DprimeVals=[DprimeVals dprime_val];
                
                abs_dprime_val = mean(abs(dprime(xinds_to_keep)));
                AbsDprimeVals=[AbsDprimeVals abs_dprime_val];
                
                DprimeStd = [DprimeStd std(dprime(xinds_to_keep))]; % std over bins, similar to abs, but here is squared error
                
                CorrCoeff = [CorrCoeff corr(ymean_hz(xinds_to_keep)', ymean_hz_base(xinds_to_keep)')]; % corr of firing rate trajcetory/
                
                % -- get median time for this bin
                median_datenum=median([segextract(trialbins).song_datenum]);
                firstday=datestr(segextract(1).song_datenum, 'ddmmmyyyy');
                eventtime=datestr(median_datenum, 'ddmmmyyyy-HHMM');
                tmp=lt_convert_EventTimes_to_RelTimes(firstday, {eventtime}); % days from start of expt
                
                MedianDatenumVals = [MedianDatenumVals median_datenum];
                MedianDayVals = [MedianDayVals tmp.FinalValue];
                
                % -- FF
                mean_FF=mean([segextract(trialbins).FF_val]);
                FFVals=[ FFVals mean_FF];
                
                % -- all trials baseline?
                datestring=WNchangeDateStrings{1};
                dnum=datenum(datestring, 'ddmmmyyyy-HHMM');
                if segextract(trialbins(end)).song_datenum < dnum
                    AreAllTrialsBaseline = [AreAllTrialsBaseline 1]; % all trials base
                    AreAllTrialsDurWN = [AreAllTrialsDurWN 0];
                elseif segextract(trialbins(1)).song_datenum >= dnum
                    AreAllTrialsBaseline = [AreAllTrialsBaseline 0];
                    AreAllTrialsDurWN = [AreAllTrialsDurWN 1];
                else
                    AreAllTrialsBaseline = [AreAllTrialsBaseline 0];
                    AreAllTrialsDurWN = [AreAllTrialsDurWN 0];
                end
            end
            
            
            % ==== PLOT
            if strcmp(typetoplot, 'mean'); % set as 1 to plot mean dprime, otherwise plots mean abs(dprime)
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots1=[hsplots1 hsplot];
                title(motif_regexpr_str{m}); ylabel('mean premotor dprime');
                plot(MedianDayVals, DprimeVals, 'ok');
                plot(MedianDayVals(logical(AreAllTrialsBaseline)), DprimeVals(logical(AreAllTrialsBaseline)), 'ob');
                plot(MedianDayVals(logical(AreAllTrialsDurWN)), DprimeVals(logical(AreAllTrialsDurWN)), 'or');
                
            elseif strcmp(typetoplot, 'abs');
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots1=[hsplots1 hsplot];
                title(motif_regexpr_str{m}); ylabel('mean abs(premotor dprime)');
                plot(MedianDayVals, AbsDprimeVals, 'ok');
                plot(MedianDayVals(logical(AreAllTrialsBaseline)), AbsDprimeVals(logical(AreAllTrialsBaseline)), 'ob');
                plot(MedianDayVals(logical(AreAllTrialsDurWN)), AbsDprimeVals(logical(AreAllTrialsDurWN)), 'or');
            elseif strcmp(typetoplot, 'std');
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots1=[hsplots1 hsplot];
                title(motif_regexpr_str{m}); ylabel('std of dprime (over time bins)');
                plot(MedianDayVals, DprimeStd, 'ok');
                plot(MedianDayVals(logical(AreAllTrialsBaseline)), DprimeStd(logical(AreAllTrialsBaseline)), 'ob');
                plot(MedianDayVals(logical(AreAllTrialsDurWN)), DprimeStd(logical(AreAllTrialsDurWN)), 'or');
            elseif strcmp(typetoplot, 'corr');
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots1=[hsplots1 hsplot];
                title(motif_regexpr_str{m}); ylabel('corr coeff of smooth fr');
                plot(MedianDayVals, CorrCoeff, 'ok');
                plot(MedianDayVals(logical(AreAllTrialsBaseline)), CorrCoeff(logical(AreAllTrialsBaseline)), 'ob');
                plot(MedianDayVals(logical(AreAllTrialsDurWN)), CorrCoeff(logical(AreAllTrialsDurWN)), 'or');                
            end
            
            if         UseEntireBaseline==0
                % -- note which bins overlap with baseline bin (i.e. first bin)
                boundaryAfterBase=mean(MedianDayVals(TrialBinSize:TrialBinSize+1));
                line([boundaryAfterBase boundaryAfterBase], ylim, 'Color','k', 'LineStyle', '--');
            end
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots2=[hsplots2 hsplot];
            title(motif_regexpr_str{m}); ylabel('ff (hz)');
            plot(MedianDayVals, FFVals, 'sk');
            plot(MedianDayVals(logical(AreAllTrialsBaseline)), FFVals(logical(AreAllTrialsBaseline)), 'sb');
            plot(MedianDayVals(logical(AreAllTrialsDurWN)), FFVals(logical(AreAllTrialsDurWN)), 'sr');
            
            
        end
        
        linkaxes(hsplots1, 'xy');
        linkaxes(hsplots2, 'x');
        lt_subtitle(['neuron' num2str(i)]);
    end
end

